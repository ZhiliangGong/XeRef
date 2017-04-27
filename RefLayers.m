classdef RefLayers < handle
    
    properties
        
        qoff = 0
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
        inhomo = 0
        
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
            
            this.qoff = paras.p0(1);
            this.ed = paras.ed;
            this.thickness = paras.thickness;
            
            if paras.pro
                this.protein = paras.protein;
                this.density = paras.density;
                this.insertion = paras.insertion;
                this.theta = paras.theta;
                this.phi = paras.phi;
                this.inhomo = paras.inhomo;
            else
                this.protein = [];
            end
            
            this.getSmoothEdProfile();
            
        end
        
        function updateRoughness(this, roughness)
            
            this.sigma = roughness;
            this.getSmoothEdProfile();
            
        end
        
        % calculation
        
        function ref = getRef(this, q)
            
            if nargin == 1
                q = this.defaultQ();
            end
            
            qc = this.getQc();
            
            ref = this.parratt(this.profile.ed, this.profile.thickness, q + this.qoff, qc);
            
            if this.inhomo > 0 && this.inhomo <= 1
%                 ed_lipid = [0.335 0.41 0.26 0];
%                 thick_lipid = [Inf 8 14.36 Inf];
                ed_lipid = [0.335 0.4427 0.3008 0];
                thick_lipid = [Inf 8 14.36 Inf];
                profile_lipid = this.calculateSmoothEdProfile(ed_lipid, thick_lipid, this.sigma, 0);
                ref_lipid = this.parratt(profile_lipid.ed, profile_lipid.thickness, q + this.qoff, qc);
                ref = this.inhomo * ref_lipid + (1 - this.inhomo) * ref;
            end
            
        end
        
        function dat = getFNR(this, q)
            
            if nargin == 1
                q = this.defaultQ();
            end
            
            dat.q = q;
            dat.ref = this.getRef(q + this.qoff);
            dat.fnr = dat.ref ./ this.getFresnel(q + this.qoff);
            
        end
        
        % fitting
        
        function fitDataQuick(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            if sum(sel) > 0
                this.fitAll(refData, para_all, lb_all, ub_all);
%                 this.fitParaOneByOne(refData, para_all, lb_all, ub_all, n_steps);
            end
            
        end
        
        function fitDataThorough(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            if sum(sel) > 0
                this.fitAll(refData, para_all, lb_all, ub_all);
                this.fitParaOneByOne(refData, para_all, lb_all, ub_all, n_steps);
%                 if sum(sel) > 1
%                     this.fitParaByPairs(refData, para_all, lb_all, ub_all, n_steps);
%                 end
            end
            
        end
        
        function fitAngleGrid(this, refData, para_all, lb_all, ub_all, n_steps)
            
            if ~ isempty(this.protein) && sum(lb_all(end-2:end-1) == ub_all(end-2:end-1)) == 0
                ind1 = length(para_all) - 2;
                ind2 = ind1 + 1;
                if isempty(this.fits)
                    this.fitAll(refData, para_all, lb_all, ub_all);
                end
                this.fits.angles = this.fitPair(refData.q, refData.ref, refData.err, this.sigma, this.energy, this.protein.pdb, this.protein.gridSize, para_all, lb_all, ub_all, ind1, ind2, n_steps);
            end
            
        end
        
        function fitAll(this, refData, para_all, lb_all, ub_all)
            
            pro_flag = ~ isempty(this.protein);
            
            if pro_flag
                n_layer = (length(para_all) - 4) / 2;
                pdb = this.protein.pdb;
                gs = this.protein.gridSize;
            else
                n_layer = (length(para_all) + 1) / 2;
                pdb = [];
                gs = [];
            end
            
            sel = lb_all ~= ub_all;
            para_partial = para_all(sel);
            lb_partial = lb_all(sel);
            ub_partial = ub_all(sel);
            
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e25, 'MaxIter', 1e6,...
                'Algorithm','trust-region-reflective', 'Display', 'iter', 'UseParallel', true);
            
            % fit all varying parameters at once
            
            sigma_temp = this.sigma;
            energy_temp = this.energy;
            
            fitAllFun = @(p) (RefLayers.calculateRefPartialPara(refData.q, p, lb_all, ub_all, sigma_temp, energy_temp, pdb, gs) - refData.ref) ./ refData.err;
            [para_fitted, chi2] = lsqnonlin(fitAllFun, para_partial, lb_partial, ub_partial, options);
            
            if isempty(this.protein)
                pro_flag = false;
            else
                pro_flag = true;
            end
            
            this.fits.all.para_names_all = this.getParaNames(n_layer, pro_flag);
            
            this.fits.all.para_names_fitted = this.fits.all.para_names_all(sel);
            para_all = lb_all;
            para_all(sel) = para_fitted;
            this.fits.all.para_all = para_all;
            this.fits.all.para_fitted = para_fitted;
            this.fits.all.fitted = sel;
            this.fits.all.chi2 = chi2;
            this.fits.all.q = refData.q;
            this.fits.all.ref = refData.ref;
            this.fits.all.err = refData.err;
            this.fits.all.ref_fit = RefLayers.calculateRefPartialPara(refData.q, para_fitted, lb_all, ub_all, sigma_temp, energy_temp, pdb, gs);
            qc = this.calculateQc(para_all(n_layer + 1), this.getWavelength());
            this.fits.all.ref_fit_fnr = this.fits.all.ref_fit ./ (this.calculateFresnel(refData.q, qc));
            
        end
        
        function fitParaOneByOne(this, refData, para_all, lb_all, ub_all, n_steps)
            
            sel = lb_all ~= ub_all;
            indices = find(sel);
            para_partial = para_all(~sel);
            
            M = length(para_partial);
            this.fits.one = cell(1, M);
            
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'Display', 'none', 'UseParallel', true);
            
            % prepare for parfor loop
            q_par = refData.q;
            ref_par = refData.ref;
            err_par = refData.err;
            sigma_par = this.sigma;
            energy_par = this.energy;
            
            if isempty(this.protein)
                pdb = [];
                gs = [];
            else
                pdb = this.protein.pdb;
                gs = this.protein.gridSize;
            end
            
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
                    
                    fitFun = @(p) (RefLayers.calculateRefPartialPara(q_par, p, lb_all_fix_one, ub_all_fix_one, sigma_par, energy_par, pdb, gs) - ref_par) ./ err_par;
                    
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
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'Display', 'none', 'UseParallel', true);
            
            switch num
                case 0
                    disp('Select at least 2 parameters to fit for pair fitting.');
                case 1
                    disp('Select at least 2 parameters to fit for pair fitting.');
                case 2
                    
                otherwise
                    results = cell(num, num);
                    tic;
                    
                    if isempty(this.protein)
                        pdb = [];
                        gs = [];
                    else
                        pdb = this.protein.pdb;
                        gs = this.protein.gridSize;
                    end
                    
                    for i = 1 : num - 1
                        for j = i + 1 : num
                            
                            para_range_1 = linspace(lb_all(indices(i)), ub_all(indices(i)), n_steps);
                            para_range_2 = linspace(lb_all(indices(j)), ub_all(indices(j)), n_steps);
                            
                            result.para_range_1 = para_range_1;
                            result.para_range_2 = para_range_2;
                            result.chi2 = zeros(n_steps, n_steps);
                            result.para_names = this.fits.all.para_names_all(indices([i, j]));
                            
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
                                    
                                    fitFun = @(p) (RefLayers.calculateRefPartialPara(q_par, p, lb_all_fix_two, ub_all_fix_two, sigma_par, energy_par, pdb, gs) - ref_par) ./ err_par;
                                    
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
        
        function result = fitPair(this, q, ref, err, sigma, energy, pdb, gs, para_all, lb_all, ub_all, ind1, ind2, n_steps)
            
            if ind1 == ind2
                error('the indices of the parameter pair must not be the same.');
            end
            
            result.chi2 = zeros(n_steps, n_steps);
            result.para_names = this.fits.all.para_names_all([ind1, ind2]);
            
            para_range_1 = linspace(lb_all(ind1), ub_all(ind1), n_steps);
            para_range_2 = linspace(lb_all(ind2), ub_all(ind2), n_steps);
            
            result.para_range_1 = para_range_1;
            result.para_range_2 = para_range_2;
            
            % prepare for parfor loop
            chi2_mat = zeros(n_steps, n_steps);
            
            options = this.options('iter');
            
            tic;
            parfor m = 1 : n_steps
                para_range_2_par = para_range_2;
                para_all_par = para_all;
                for n = 1 : n_steps
                    lb_all_fix_two = lb_all;
                    lb_all_fix_two(ind1) = para_range_1(m);
                    lb_all_fix_two(ind2) = para_range_2_par(n);
                    ub_all_fix_two = ub_all;
                    ub_all_fix_two(ind1) = para_range_1(m);
                    ub_all_fix_two(ind2) = para_range_2_par(n);
                    
                    fitFun = @(p) (RefLayers.calculateRefPartialPara(q, p, lb_all_fix_two, ub_all_fix_two, sigma, energy, pdb, gs) - ref) ./ err;
                    
                    newsel = lb_all_fix_two == ub_all_fix_two;
                    lb_partial_fix_two = lb_all_fix_two(~newsel);
                    ub_partial_fix_two = ub_all_fix_two(~newsel);
                    para_partial_fix_two = para_all_par(~newsel);
                    [~, chi2] = lsqnonlin(fitFun, para_partial_fix_two, lb_partial_fix_two, ub_partial_fix_two, options);
                    
                    chi2_mat(m, n) = chi2;
                    
                end
            end
            toc;
            
            lk = chi2_mat;
            lk = exp(-(lk-min(lk(:)))/2);
            result.likelihood = lk/sum(lk(:));
            result.chi2 = chi2_mat;
            result.confidence = this.confidenceContour(para_range_1, para_range_2, result.likelihood, 0.95);
            
        end
        
        function R = getFresnel(this, q)
            
            if nargin == 1
                q = this.getQ();
            end
            
            qc = this.getQc();
            
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
                [pro_ed, pro_thick, pro_area] = this.protein.getEdProfile(this.theta, this.phi);
                [com_ed, com_thick, z_origin] = this.getCompositeEdProfile(this.ed, this.thickness, pro_ed, pro_thick, pro_area, this.density, this.insertion);
                this.profile = this.calculateSmoothEdProfile(com_ed, com_thick, this.sigma, z_origin);
            end
            
        end
        
        % plot
        
        function plotEdProfile(this, ax)
            
            if nargin == 1 || isempty(ax)
                figure;
                ax = gca;
            end
            
            plot(ax, this.profile.z, this.profile.layerEd, '-k', 'linewidth', 2);
            hold(ax, 'on');
            plot(ax, this.profile.z, this.profile.ed, '-b', 'linewidth', 2);
            legend(ax, {'Layer Structure', 'Smoothed With Roughness'});
            xlabel(ax, '$$ z (\AA) $$', 'interpreter', 'latex', 'fontsize', 14);
            ylabel(ax, '$$ Electron Density (\AA^{-3}) $$', 'interpreter', 'latex', 'fontsize', 14);
            hold(ax, 'off');
            
        end
        
    end
    
    methods(Static)
        
        % ed profiles
        
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
        
        function [ed, thick, z_origin] = getCompositeEdProfile(layer_ed, layer_thick, pro_ed, pro_thick, pro_area, density, insertion)
            
            gs = pro_thick(1);
            pro_top = insertion;
            
            % top and bottom position of protein
            pro_length = sum(pro_thick);
            grid_n = length(pro_ed);
            pro_bot = pro_top - pro_length;
            
            interface_z = cumsum([0, layer_thick( 2 : end - 1)]);
            layer_top_z = [interface_z, Inf];
            layer_bot_z = [-Inf, interface_z];
            
            % obtain the layers without intersection with the protein
            sel_top = layer_bot_z >= pro_top | layer_top_z == Inf;
            ed_top = layer_ed(sel_top);
            thick_top = layer_thick(sel_top);
            
            sel_bot = layer_top_z <= pro_bot | layer_bot_z == -Inf;
            ed_bot = layer_ed(sel_bot);
            thick_bot = layer_thick(sel_bot);
            
            % find the layer that has the top/bottom edges of the protein
            sel_pro_top = layer_top_z >= pro_top & layer_bot_z < pro_top;
            sel_pro_bot = layer_top_z >= pro_bot & layer_bot_z < pro_bot;
            
            % obtain the top part of protein free part of the layer
            % that contains the top of the protein, or the other way
            thick_partial_top = layer_top_z(sel_pro_top) - pro_top;
            if thick_partial_top ~= 0 && thick_partial_top ~= Inf
                ed_top = [layer_ed(sel_pro_top), ed_top];
                thick_top = [thick_partial_top, thick_top];
            end
            
            thick_partial_bot = pro_bot - layer_bot_z(sel_pro_bot);
            if thick_partial_bot ~=0 && thick_partial_bot ~= Inf
                ed_bot = [ed_bot, layer_ed(sel_pro_bot)];
                thick_bot = [thick_bot, thick_partial_bot];
            end
            
            % find the interfaces that goes through the protein
            through_z = interface_z(interface_z < pro_top & interface_z > pro_bot);
            pro_n = length(pro_ed);
            pro_grid_top_z = (1 : pro_n) * gs - pro_length + insertion;
            pro_grid_bot_z = pro_grid_top_z - gs;
            
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
            pro_grid_top_z_new = cumsum(pro_thick_new) - pro_length + insertion;
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
                area_frac = pro_area_new(sel_between) * density / 1e4;
                composite_ed(sel_between) = pro_ed_new(sel_between) .* area_frac + layer_ed(indices(i)) * (1 - area_frac);
            end
            
            ed = [ed_bot, composite_ed, ed_top];
            thick = [thick_bot, pro_thick_new, thick_top];
            z_origin = pro_length - insertion;
            
        end
        
        % reflectivity calculation
        
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
        
        function f = calculateFresnel(q, qc)
            
            r = (q - sqrt(q.^2 - ones(size(q)) * qc^2)) ./ (q + sqrt(q.^2 - ones(size(q)) * qc^2));
            f = r .* conj(r);
            
        end
        
        function ref = calculateRefPartialPara(q, para_partial, lb_all, ub_all, sigma, energy, pdb, gs)
            
            para_full = lb_all;
            sel = lb_all == ub_all;
            para_full(~sel) = para_partial;
            
            qoff = para_full(1);
            
            pro_inhomo = 0;
            
            if isempty(pdb)
                n_layer = (length(para_full) + 1) / 2;
                thick_coarse = [Inf, para_full(2 * n_layer - 1 : -1 : n_layer + 2), Inf];
                ed_coarse = para_full(n_layer + 1 : -1 : 2);
                z_origin = 0;
            else
                n_layer = (length(para_full) - 4) / 2;
                
                layer_thick = [Inf, para_full(2 * n_layer - 1 : -1 : n_layer + 2), Inf];
                layer_ed = para_full(n_layer + 1 : -1 : 2);
                
                pro_inhomo = para_full(end);
                pro_phi = para_full(end - 1);
                pro_theta = para_full(end - 2);
                pro_insertion = para_full(end - 3);
                pro_density = para_full(end - 4);
                [pro_ed, pro_thick, pro_area] = RefProtein.calculateEdProfile(pdb.x, pdb.y, pdb.z, pdb.radius, pdb.electron, pro_theta, pro_phi, gs);
                
                [ed_coarse, thick_coarse, z_origin] = RefLayers.getCompositeEdProfile(layer_ed, layer_thick, pro_ed, pro_thick, pro_area, pro_density, pro_insertion);
                
            end
               
            prof = RefLayers.calculateSmoothEdProfile(ed_coarse, thick_coarse, sigma, z_origin);
            wl = RefLayers.calculateWavelength(energy);
            qc = RefLayers.calculateQc(ed_coarse(1), wl);
            
            ref = RefLayers.parratt(prof.ed, prof.thickness, q + qoff, qc);
            
            if pro_inhomo > 0 && pro_inhomo <= 1
                ed_lipid = [0.335 0.4427 0.3008 0];
                thick_lipid = [Inf 8 14.36 Inf];
                profile_lipid = RefLayers.calculateSmoothEdProfile(ed_lipid, thick_lipid, sigma, 0);
                ref_lipid = RefLayers.parratt(profile_lipid.ed, profile_lipid.thickness, q + qoff, qc);
                ref = pro_inhomo * ref_lipid + (1 - pro_inhomo) * ref;
            elseif pro_inhomo > 1
                warning('The inhomogeneity parameter must be between 0 and 1. Ignoring this parameter.');
            end
            
        end
        
        % utility
        
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
        
        function names = getParaNames(n_layer, pro_flag)
            
            if pro_flag
                names = cell(1, n_layer * 2 + 4);
                names(end -4 : end) = {'Density'; 'Insertion'; 'Theta'; 'Phi'; 'Inhomogeneity'};
%                 names{end} = 'Inhomogeneity';
%                 names{end - 1} = 'Phi';
%                 names{end - 2} = 'Theta';
%                 names{end - 3} = 'Insertion';
%                 names{end - 4} = 'Density';
            else
                names = cell(1, n_layer * 2 - 1);
            end
            
            names{1} = 'Qz-Offset';
            names{2} = 'Top-ED';
            names{n_layer + 1} = 'Bottom-ED';
            for i = 1 : n_layer - 2
                names{n_layer + 1 - i} = ['Layer-', num2str(i), '-ED'];
                names{2 * n_layer - i} = ['Layer-', num2str(i), '-Thkns'];
            end
            
        end
        
        function q = defaultQ()
            
            q = 0.03 : 0.005 : 0.7;
            
        end
        
        function stepSizes = getStepSizes(n_layer, pro_flag)
            
            if pro_flag
                stepSizes = [1e-6, ones(1, n_layer) * 1e-4, ones(1, n_layer - 2) * 1e-2, 1e-2, 1e-2, 1, 1];
            else
                stepSizes = [1e-6, ones(1, n_layer) * 1e-4, ones(1, n_layer - 2) * 1e-3];
            end
            
        end
        
        % constants
        
        function o = options(display)
            
            if nargin == 0
                display = 'off';
            end
            
            switch display
                case 'iter'
                    o = optimoptions('lsqnonlin', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'Display', 'iter', 'UseParallel', true);
                case 'off'
                    o = optimoptions('lsqnonlin', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'Display', 'none', 'UseParallel', true);
            end
            
        end
        
    end
    
end