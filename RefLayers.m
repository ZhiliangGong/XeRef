classdef RefLayers < handle
    
    properties
        
        ed
        thickness
        
        fits
        
        profile
        sigma = 3.3
        energy
        
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
            
            this.profile = this.generateSmoothProfile(this.ed, this.thickness, this.sigma);
            
        end
        
        % calculation
        
        function ref = getRef(this, q)
            
            if nargin == 1
                q = this.defaultQ();
            end
            
            qc = this.getQc();
            
            ref = this.parratt(this.profile.ed, this.profile.thickness, q, qc);
            
        end
        
        function ref = getRefForFullPara(this, q, paras)
            
            n_layer = (length(paras) - 3) / 2;
            qoff = paras(1);
            thick = [Inf, paras(2 : n_layer + 1), Inf];
            ED = paras(n_layer + 2 : end);
            
            pro = this.generateSmoothProfile(ED, thick, this.sigma);
            qc = this.calculateQc(ED(1), this.getWavelength());
            
            ref = this.parratt(pro.ed, pro.thickness, q + qoff, qc);
            
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
        
        function fitData(this, refData, para_all, lb_all, ub_all, steps)
            
            n_layer = (length(para_all) + 1) / 2;
            sel = lb_all == ub_all;
            para_partial = para_all(~sel);
            lb_partial = lb_all(~sel);
            ub_partial = ub_all(~sel);
            
            options = optimoptions('lsqnonlin', 'MaxFunEvals', 1e25, 'MaxIter', 1e5, 'Display', 'iter');
            fitFun = @(p) (this.getRefForPartialPara(refData.q, p, lb_all, ub_all) - refData.ref) ./ refData.err;
            [para_fitted, chi2] = lsqnonlin(fitFun, para_partial, lb_partial, ub_partial, options);
            
            this.fits.all.para_names_all = this.generateParaNames(n_layer);
            this.fits.all.para_names_fitted = this.fits.all.para_names_all(~sel);
            para_all = lb_all;
            para_all(~sel) = para_fitted;
            this.fits.all.para_all = para_all;
            this.fits.all.para_fitted = para_fitted;
            this.fits.all.chi2 = chi2;
            this.fits.all.q = refData.q;
            this.fits.all.ref = refData.ref;
            this.fits.all.err = refData.err;
            this.fits.all.ref_fit = this.getRefForPartialPara(refData.q, para_fitted, lb_all, ub_all);
            qc = this.calculateQc(para_all(end - n_layer + 1), this.getWavelength());
            this.fits.all.ref_fit_fnr = this.fits.all.ref_fit ./ (this.calculateFresnel(refData.q + para_all(1), qc));
            
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
        
        function profile = getSmoothEdProfile(this, sigma)
            
            if nargin == 2
                this.sigma = sigma;
            end
            
            profile = this.generateSmoothProfile(this.ed, this.thickness, this.sigma);
            
        end
        
    end
    
    methods(Static)
        
        function profile = generateSmoothProfile(ed, thickness, sigma)
            
            if nargin == 2
                sigma = 3.3;
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
            profile.z = binPos';
            profile.thickness = ones(1, m) * binSize;
            ed_mat = repmat(ed, m);
            profile.layerEd = ed_mat(indices)';
            
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
            names{n_layer} = 'Bottom-ED';
            names{n_layer * 2 - 1} = 'Top_ED';
            for i = 1 : n_layer - 2
                names{i + 1} = ['Layer-', num2str(i), '-Thkns'];
                names{n_layer + i} = ['Layer-', num2str(i), '-ED'];
            end
            
        end
        
        function f = calculateFresnel(q, qc)
            
            r = (q - sqrt(q.^2 - ones(size(q)) * qc^2)) ./ (q + sqrt(q.^2 - ones(size(q)) * qc^2));
            f = r .* conj(r);
            
        end
        
    end
    
end