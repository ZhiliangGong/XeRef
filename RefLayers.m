classdef RefLayers < handle
    
    properties
        
        ed
        thickness
        
        fits
        
        profile
        sigma = 3.4
        energy
        
        protein
        
        density = 0
        insertion = 0
        theta = 0
        phi = 0
        
    end
    
    methods
        
        function this = RefLayers(energy, ed, thickness, sigma)
            
            this.energy = energy;
            if nargin >= 3
                this.ed = ed;
                this.thickness = thickness;
                if nargin == 4
                    this.sigma = sigma;
                end
            end
            
            this.getSmoothEdProfile();
            
        end
        
        function updateModel(this, paras)
            
            this.ed = paras.ed;
            this.thickness = paras.thickness;
            
            if paras.pro
                this.protein = paras.protein;
                this.density = paras.density;
                this.insertion = paras.insertion;
                this.theta = paras.theta;
                this.phi = paras.phi;
            else
                this.protein = [];
            end
            
            this.getSmoothEdProfile();
            
        end
        
        % calculation
        
        function ref = getRef(this, q)
            
            if nargin == 1
                q = this.defaultQ();
            end
            
            qc = this.getQc();
            
            ref = this.parratt(this.profile.ed, this.profile.thickness, q, qc);
            
        end
        
        function ref = getRefForFullPara(this, q, full_para)
            
            if isempty(this.protein)
                n_layer = (length(full_para) + 1) / 2;
            else
                n_layer = (length(full_para) - 3) / 2;
            end
            qoff = full_para(1);
            thick = [Inf, full_para(n_layer * 2 - 1: -1 : n_layer + 2), Inf];
            ED = full_para(n_layer + 1 : -1 : 2);
            
            edProfile = this.calculateSmoothEdProfile(ED, thick, this.sigma);
            qc = this.calculateQc(ED(1), this.getWavelength());
            
            ref = this.parratt(edProfile.ed, edProfile.thickness, q + qoff, qc);
            
        end
        
        function ref = getRefForPartialPara(this, q, para_partial, lb_all, ub_all)
            
            para_full = lb_all;
            sel = lb_all == ub_all;
            para_full(~sel) = para_partial;
            
            ref = this.getRefForFullPara(q, para_full);
            
        end
        
        function dat = getFNR(this, q)
            
            if nargin == 1
                q = this.defaultQ();
            end
            
            dat.q = q;
            dat.ref = this.getRef(q);
            dat.fnr = dat.ref ./ this.getFresnel(q);
            
        end
        
        % fitting
        
        function fitDataQuick(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            if sum(sel) > 0
                this.fitAll(refData, para_all, lb_all, ub_all);
                this.fitParaOneByOne(refData, para_all, lb_all, ub_all, n_steps);
            end
            
        end
        
        function fitDataThorough(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            if sum(sel) > 0
                this.fitAll(refData, para_all, lb_all, ub_all);
                this.fitParaOneByOne(refData, para_all, lb_all, ub_all, n_steps);
                if sum(sel) > 1
                    this.fitParaByPairs(refData, para_all, lb_all, ub_all, n_steps);
                end
            end
            
        end
        
        function fitAll(this, refData, para_all, lb_all, ub_all)
            
            n_layer = (length(para_all) + 1) / 2;
            sel = lb_all == ub_all;
            para_partial = para_all(~sel);
            lb_partial = lb_all(~sel);
            ub_partial = ub_all(~sel);
            
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e25, 'MaxIter', 1e5, 'Display', 'none');
            
            % fit all varying parameters at once
            
            fitAllFun = @(p) (this.getRefForPartialPara(refData.q, p, lb_all, ub_all) - refData.ref) ./ refData.err;
            [para_fitted, chi2] = lsqnonlin(fitAllFun, para_partial, lb_partial, ub_partial, options);
            
            this.fits.all.para_names_all = this.generateParaNames(n_layer);
            this.fits.all.para_names_fitted = this.fits.all.para_names_all(~sel);
            para_all = lb_all;
            para_all(~sel) = para_fitted;
            this.fits.all.para_all = para_all;
            this.fits.all.para_fitted = para_fitted;
            this.fits.all.fitted = ~sel;
            this.fits.all.chi2 = chi2;
            this.fits.all.q = refData.q;
            this.fits.all.ref = refData.ref;
            this.fits.all.err = refData.err;
            this.fits.all.ref_fit = this.getRefForPartialPara(refData.q, para_fitted, lb_all, ub_all);
            qc = this.calculateQc(para_all(n_layer + 1), this.getWavelength());
            this.fits.all.ref_fit_fnr = this.fits.all.ref_fit ./ (this.calculateFresnel(refData.q, qc));
            
        end
        
        function fitParaOneByOne(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            indices = find(sel);
            para_partial = para_all(~sel);
            
            M = length(para_partial);
            this.fits.one = cell(1, M);
            
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e25, 'MaxIter', 1e5, 'Display', 'none');
            
            % prepare for parfor loop
            q_par = refData.q;
            ref_par = refData.ref;
            err_par = refData.err;
            sigma_par = this.sigma;
            energy_par = this.energy;
            
            tic;
            for i = 1 : length(indices)
                ind = indices(i);
                para_range = linspace(lb_all(ind), ub_all(ind), n_steps);
                    
                result.para_range = para_range;
                chi2_array = zeros(1, n_steps);
                result.para_name = this.fits.all.para_names_all{ind};
                
                parfor j = 1 : n_steps
                    
                    para_all_par = para_all;
                    
                    lb_all_fix_one = lb_all;
                    lb_all_fix_one(ind) = para_range(j);
                    ub_all_fix_one = ub_all;
                    ub_all_fix_one(ind) = para_range(j);
                    
                    fitFun = @(p) (RefLayers.calculateRefPartialPara(q_par, p, lb_all_fix_one, ub_all_fix_one, sigma_par, energy_par) - ref_par) ./ err_par;
                    
                    newsel = lb_all_fix_one ~= ub_all_fix_one;
                    lb_partial_fix_one = lb_all_fix_one(newsel);
                    ub_partial_fix_one = ub_all_fix_one(newsel);
                    para_partial_fix_one = para_all_par(newsel);
                    
                    [~, chi2] = lsqnonlin(fitFun, para_partial_fix_one, lb_partial_fix_one, ub_partial_fix_one, options);
                    chi2_array(j) = chi2;
                    
                end
                
                result.chi2 = chi2_array;
                lk = exp( - (result.chi2 - min(result.chi2)) / 2);
                result.likelihood = lk / sum(lk);
                
                [gauss_para, flag] = fitGaussianLikelihood(result.para_range, result.likelihood);
                if flag
                    warning('%s %s', result.para_name, 'fitting bad.');
                end
                result.gauss.x = linspace(para_range(1), para_range(end), 100);
                result.gauss.y = this.bellCurve(result.gauss.x, gauss_para);
                result.gauss.para = gauss_para;
                
                this.fits.one{i} = result;
            end
            toc;
            
        end
        
        function fitParaByPairs(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            indices = find(sel);
            num = length(indices);
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e25, 'MaxIter', 1e5, 'Display', 'none');
            
            switch num
                case 0
                    disp('Select at least 2 parameters to fit for pair fitting.');
                case 1
                    disp('Select at least 2 parameters to fit for pair fitting.');
                case 2
                    
                otherwise
                    results = cell(num, num);
                    tic;
                    for i = 1 : num - 1
                        for j = i + 1 : num
                            para_range_1 = linspace(lb_all(indices(i)), ub_all(indices(i)), n_steps);
                            para_range_2 = linspace(lb_all(indices(j)), ub_all(indices(j)), n_steps);
                            
                            result.para_range_1 = para_range_1;
                            result.para_range_2 = para_range_2;
                            result.chi2 = zeros(n_steps, n_steps);
                            result.para_names = this.fits.all.para_names_all([i, j]);
                            
                            % prepare for parfor loop
                            q_par = refData.q;
                            ref_par = refData.ref;
                            err_par = refData.err;
                            sigma_par = this.sigma;
                            energy_par = this.energy;
                            chi2_mat = zeros(n_steps, n_steps);
                            
                            parfor m = 1 : n_steps
                                para_range_2_par = para_range_2;
                                para_all_par = para_all;
                                indices_par = indices;
                                for n = 1 : n_steps
                                    lb_all_fix_two = lb_all;
                                    lb_all_fix_two(indices_par(i)) = para_range_1(m);
                                    lb_all_fix_two(indices_par(j)) = para_range_2_par(n);
                                    ub_all_fix_two = ub_all;
                                    ub_all_fix_two(indices_par(i)) = para_range_1(m);
                                    ub_all_fix_two(indices_par(j)) = para_range_2_par(n);
                                    
                                    fitFun = @(p) (RefLayers.calculateRefPartialPara(q_par, p, lb_all_fix_two, ub_all_fix_two, sigma_par, energy_par) - ref_par) ./ err_par;
                                    
                                    newsel = lb_all_fix_two == ub_all_fix_two;
                                    lb_partial_fix_two = lb_all_fix_two(~newsel);
                                    ub_partial_fix_two = ub_all_fix_two(~newsel);
                                    para_partial_fix_two = para_all_par(~newsel);
                                    [~, chi2] = lsqnonlin(fitFun, para_partial_fix_two, lb_partial_fix_two, ub_partial_fix_two, options);
                                    
                                    chi2_mat(m, n) = chi2;
                                    
                                end
                            end
                            
                            lk = chi2_mat;
                            lk = exp(-(lk-min(lk(:)))/2);
                            result.likelihood = lk/sum(lk(:));
                            result.chi2 = chi2_mat;
                            result.confidence = this.confidenceContour(para_range_1, para_range_2, result.likelihood, 0.95);
                            
                            results{i, j} = result;
                            
                        end
                    end
                    toc;
                    this.fits.two = results;
            end
            
        end
        
        % utility
        
        function R = getFresnel(this, q)
            
            if nargin == 1
                q = this.getQ();
            end
            
            qc = this.getQc();
            
%             r = (q - sqrt(q.^2 - ones(size(q)) * qc^2)) ./ (q + sqrt(q.^2 - ones(size(q)) * qc^2));
%             R = r .* conj(r);
            
            R = this.calculateFresnel(q, qc);
            
        end
        
        function qc = getQc(this)
            
            qc = this.calculateQc(this.ed(1), this.getWavelength());
            
        end
        
        function wl = getWavelength(this)
            
            wl = this.calculateWavelength(this.energy);
            
        end
        
        function getSmoothEdProfile(this)
            
            if isempty(this.protein) || this.density == 0
                this.profile = this.calculateSmoothEdProfile(this.ed, this.thickness, this.sigma, 0);
            else
                [ed_pro, thick_pro, z_origin] = this.getCompositeEdProfile();
                this.profile = this.calculateSmoothEdProfile(ed_pro, thick_pro, this.sigma, z_origin);
            end
            
        end
        
    end
    
    methods(Static)
        
        function profile = calculateSmoothEdProfile(ed, thickness, sigma, z_origin)
            
            if isempty(z_origin)
                z_origin = 0;
            end
            
            binSize = 0.25;
            transSigma = 7;
            transThickness = sum(thickness(2 : end - 1)) + 2 * transSigma * sigma;
            binNumber = ceil(transThickness / binSize);
            
            thick = [0, thickness(2:end-1), 0];
            
            layerCenterPos = cumsum(thick) - thick / 2;
            binPos = cumsum(ones(binNumber, 1) * binSize) - (transSigma * sigma) - binSize / 2;
            
            m = length(binPos);
            n = length(layerCenterPos);
            contributions = zeros(m, n);
            
            if n > 2
                binPos_mat = repmat(binPos, 1, n - 2);
                layerCenterPos_mat = repmat(layerCenterPos(2:end-1), m, 1);
                thick_mat = repmat(thick(2:end-1), m, 1);
                ed_mat = repmat(ed(2:end-1), m, 1);
                contributions(:, 2: end-1) = (erf((binPos_mat - layerCenterPos_mat + thick_mat / 2) / sqrt(2) / sigma) ...
                    + erf((- binPos_mat + layerCenterPos_mat + thick_mat / 2) / sqrt(2) / sigma)) / 2 .* ed_mat ;
            end
            
            contributions(:, 1) = ed(1) * (erf(-binPos / sqrt(2) / sigma) + 1) / 2;
            
            breakPoints = [-Inf, 0, cumsum(thickness(2:end))];
            
            indices = binPos < repmat(breakPoints(2:end), m, 1) & binPos > repmat(breakPoints(1:end-1), m, 1);
            
            binEd = sum(contributions, 2);
            profile.ed = binEd';
            profile.z = binPos' - z_origin;
            profile.thickness = ones(1, m) * binSize;
            ed_mat = repmat(ed, m);
            profile.layerEd = ed_mat(indices)';
            
        end
        
        function [ed, thick, z_origin] = getCompositeEdProfile(layer_ed, layer_thick, pro_ed, pro_thick, pro_area)
            
            if isempty(this.protein)
                ed = this.ed;
                thick = this.thickness;
                z_origin = 0;
            else
                [pro_ed, pro_thick, pro_area] = this.protein.generateSingleEdProfile(this.theta, this.phi);
                pro_top = this.insertion;
                num_per_100nm2 = this.density;
                
                % top and bottom position of protein
                pro_length = sum(pro_thick);
                grid_n = length(pro_ed);
                pro_bot = pro_top - pro_length;
                
                interface_z = cumsum([0, this.thickness( 2 : end - 1)]);
                layer_top_z = [interface_z, Inf];
                layer_bot_z = [-Inf, interface_z];
                
                % obtain the layers without intersection with the protein
                sel_top = layer_bot_z >= pro_top | layer_top_z == Inf;
                ed_top = this.ed(sel_top);
                thick_top = this.thickness(sel_top);
                
                sel_bot = layer_top_z <= pro_bot | layer_bot_z == -Inf;
                ed_bot = this.ed(sel_bot);
                thick_bot = this.thickness(sel_bot);
                
                % find the layer that has the top/bottom edges of the protein
                sel_pro_top = layer_top_z >= pro_top & layer_bot_z < pro_top;
                sel_pro_bot = layer_top_z >= pro_bot & layer_bot_z < pro_bot;
                
                % obtain the top part of protein free part of the layer
                % that contains the top of the protein, or the other way
                thick_partial_top = layer_top_z(sel_pro_top) - pro_top;
                if thick_partial_top ~= 0 && thick_partial_top ~= Inf
                    ed_top = [this.ed(sel_pro_top), ed_top];
                    thick_top = [thick_partial_top, thick_top];
                end
                
                thick_partial_bot = pro_bot - layer_bot_z(sel_pro_bot);
                if thick_partial_bot ~=0 && thick_partial_bot ~= Inf
                    ed_bot = [ed_bot, this.ed(sel_pro_bot)];
                    thick_bot = [thick_bot, thick_partial_bot];
                end
                
                % find the interfaces that goes through the protein
                through_z = interface_z(interface_z < pro_top & interface_z > pro_bot);
                pro_n = length(pro_ed);
                pro_grid_top_z = (1 : pro_n) * this.protein.gridSize - pro_length + this.insertion;
                pro_grid_bot_z = pro_grid_top_z - this.protein.gridSize;
                
                % find the grids to be broken into two
                sel_break = through_z' <= pro_grid_top_z & through_z' > pro_grid_bot_z;
                [~, ind_2] = find(sel_break);
                
                sel_ed = sort([1 : grid_n, ind_2']);
                pro_ed_new = pro_ed(sel_ed);
                pro_area_new = pro_area(sel_ed);
                pro_thick_new = pro_thick(sel_ed);
                
                for i = 1 : length(ind_2)
                    pro_thick_new(ind_2(i) + i - 1) = through_z(i) - pro_grid_bot_z(ind_2(i));
                    pro_thick_new(ind_2(i) + i ) = pro_grid_top_z(ind_2(i)) - through_z(i);
                end
                
                sel_nonzero = pro_thick_new ~= 0;
                pro_thick_new = pro_thick_new(sel_nonzero);
                pro_ed_new = pro_ed_new(sel_nonzero);
                pro_area_new = pro_area_new(sel_nonzero);
                pro_grid_top_z_new = cumsum(pro_thick_new) - pro_length + this.insertion;
                pro_grid_bot_z_new = pro_grid_top_z_new - pro_thick_new;
                
                % identify the layers that intersets with the protein and
                % loop through them
                sel_overlap = (layer_top_z >= pro_top & layer_bot_z < pro_top)...
                    | (layer_top_z >= pro_bot & layer_bot_z < pro_bot)...
                    | (layer_top_z < pro_top & layer_bot_z > pro_bot);
                indices = find(sel_overlap);
                composite_ed = zeros(size(pro_ed_new));
                for i = 1 : length(indices)
                    sel_between = pro_grid_bot_z_new >= layer_bot_z(indices(i)) & pro_grid_top_z_new <= layer_top_z(indices(i));
                    area_frac = pro_area_new(sel_between) * num_per_100nm2 / 1e4;
                    composite_ed(sel_between) = pro_ed_new(sel_between) .* area_frac + this.ed(indices(i)) * (1 - area_frac);
                end
                
                ed = [ed_bot, composite_ed, ed_top];
                thick = [thick_bot, pro_thick_new, thick_top];
                z_origin = pro_length - this.insertion;
                
            end
            
        end
        
        function ref = parratt(ed, thickness, q, qc)
            
            normEd = ed / ed(1);
            thick = [0, thickness(2:end-1), 0];
            
            r = (sqrt(q.^2 - ones(size(q)) * normEd(2) * qc^2) - sqrt(q.^2 - ones(size(q)) * normEd(1) * qc^2)) ./...
                (sqrt(q.^2 - ones(size(q)) * normEd(2) * qc^2) + sqrt(q.^2 - ones(size(q)) * normEd(1) * qc^2));
            
            qj = sqrt(q.^2 - ones(size(q)) * normEd(2) * qc^2);
            
            for j = 2 : (length(normEd) - 1)
                
                qjp1 = sqrt(q.^2 - ones(size(q)) * normEd(j+1) * qc^2);
                reff = (qjp1 - qj)./(qjp1 + qj);
                phase = exp(1i * qj * thick(j));
                n1 = r .* phase;
                r = (reff + n1) ./ (1 + reff .* n1);
                
                qj = qjp1;
                
            end
            
            ref = r .* conj(r);
            
        end
        
        function q = defaultQ()
            
            q = 0.026 : 0.005 : 0.7;
            
        end
        
        function wl = calculateWavelength(energy)
            
            c = 299792458;
            h = 6.62607004e-34;
            ev2j = 1.60218e-19;
            
            wl = h * c / (energy * 1000 * ev2j) * 1e10;
            
        end
        
        function qc = calculateQc(bufferEd, wavelength)
            
            ro = 2.818*10^-5;
            k = 2 * pi / wavelength;
            delta = 2 * pi * ro * bufferEd / k^2;
            qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        function names = generateParaNames(n_layer)
            
            names = cell(1, n_layer * 2 - 1);
            names{1} = 'Qz-Offset';
            names{2} = 'Top-ED';
            names{n_layer + 1} = 'Bottom-ED';
            for i = 1 : n_layer - 2
                names{n_layer + 1 - i} = ['Layer-', num2str(i), '-ED'];
                names{2 * n_layer - i} = ['Layer-', num2str(i), '-Thkns'];
            end
            
        end
        
        function f = calculateFresnel(q, qc)
            
            r = (q - sqrt(q.^2 - ones(size(q)) * qc^2)) ./ (q + sqrt(q.^2 - ones(size(q)) * qc^2));
            f = r .* conj(r);
            
        end
        
        function ref = calculateRefPartialPara(q, para_partial, lb_all, ub_all, sigma, energy)
            
            para_full = lb_all;
            sel = lb_all == ub_all;
            para_full(~sel) = para_partial;
            
            n_layer = (length(para_full) + 1) / 2;
            qoff = para_full(1);
            thick = [Inf, para_full(end : -1 : end - n_layer + 3), Inf];
            ED = fliplr(para_full(2 : n_layer + 1));
            
            pro = RefLayers.calculateSmoothEdProfile(ED, thick, sigma);
            wl = RefLayers.calculateWavelength(energy);
            qc = RefLayers.calculateQc(ED(1), wl);
            
            ref = RefLayers.parratt(pro.ed, pro.thickness, q + qoff, qc);
            
        end
        
        function y = bellCurve(x, paras)
            
            y = paras(1) * exp( -(x - paras(2)).^2 / 2 / paras(3)^2 );
            
        end
        
        function result = confidenceContour(xdata, ydata, likelihood, confidence)
            
            if nargin == 3
                confidence = 0.95;
            end
            
            lk1 = sort(likelihood(:));
            lksum = cumsum(lk1);
            lksum = abs(lksum - (1 - confidence));
            ind = find(lksum == min(lksum), 1);
            cLevel = lk1(ind);
            
            C = contourc(xdata,ydata,likelihood,[cLevel,cLevel]);
            C = C(:,C(1,:) >= min(xdata));
            C = C(:,C(1,:) <= max(xdata));
            C = C(:,C(2,:) >= min(ydata));
            C = C(:,C(2,:) <= max(ydata));
            
            [ind1, ind2] = find(likelihood == max(likelihood(:)), 1);
            centerIndices = [ind1, ind2];
            centerValues = [xdata(ind1), ydata(ind2)];
            confidenceWindow = [min(C(1,:)),max(C(1,:));min(C(2,:)),max(C(2,:))];
            
            result.contour = C;
            result.centerIndices = centerIndices;
            result.centerValues = centerValues;
            result.confidenceWindow = confidenceWindow;
            
        end
        
    end
    
end