classdef RefLayers < handle
    
    properties
        
        ed = [0, 0.335]
        thickness = [Inf, Inf]
        
        qoff = 0
        data
        energy
        
    end
    
    methods
        
        function this = RefLayers(rawdata)
            
            this.energy = rawdata.energy;
            this.data = rawdata.goodData();
            this.qoffFit();
            
        end
        
        function qoffFit(this)
            
            qOffFitFunction = @(qoff) (this.getRefWithQoff(qoff) - this.data.ref) ./ this.data.err;
            this.qoff = lsqnonlin(qOffFitFunction, 0);
            
        end
        
        function ref = getRefWithQoff(this, qoff)
            
            ref = this.getRef(this.data.q + qoff);
            
        end
        
        function ref = getRef(this, q)
            
            if nargin == 1
                q = this.getQ();
            end
            
            qc = this.getQc();
            normEd = this.ed / this.ed(end);
            
            r = (sqrt(q.^2 - ones(size(q)) * normEd(end-1) * qc^2) - sqrt(q.^2 - ones(size(q)) * normEd(end) * qc^2)) ./...
                (sqrt(q.^2 - ones(size(q)) * normEd(end-1) * qc^2) + sqrt(q.^2 - ones(size(q)) * normEd(end) * qc^2));
            
            qj = sqrt(q.^2 - ones(size(q)) * normEd(end-1) * qc^2);
            
            for j = (length(normEd) - 1) : -1 : 2
                
                qjp1 = sqrt(q.^2 - ones(size(q)) * normEd(j-1) * qc^2);
                reff = (qjp1 - qj)./(qjp1 + qj);
                phase = exp(1i * qj * ddlay(j));
                n1 = r .* phase;
                r = (reff + n1)./(1 + reff .* n1);
                
                if abs(normEd(j - 1) - normEd(1)) < 1e-7
                    break;
                end
                
                qj = qjp1;
                
            end
            
            ref = r .* conj(r);
            
        end
        
        function fnr = getFNR(this, q)
            
            if nargin == 1
                q = this.getQ();
            end
            
            fnr = this.getRef(q) ./ this.getFresnel(q);
            
        end
        
        function d = getFresnelNormalizedData(this)
            
            fresnel = this.getFresnel();
            d.ref = this.data.ref ./ fresnel;
            d.err = this.data.err ./ fresnel;
            
        end
        
        % utility
        
        function R = getFresnel(this, q)
            
            if nargin == 1
                q = this.getQ();
            end
            
            qc = this.getQc();
            
            r = (q - sqrt(q.^2 - ones(size(q)) * qc^2)) ./ (q + sqrt(q.^2 - ones(size(q)) * qc^2));
            R = r .* conj(r);
            
        end
        
        function qc = getQc(this)
            
            ro = 2.818*10^-5; % radius of electron in angstrom
            k = 2 * pi / this.getWavelength();
            delta = 2 * pi * ro * this.ed(end) / k^2;
            qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        function wl = getWavelength(this)
            
            c = 299792458;
            h = 6.62607004e-34;
            ev2j = 1.60218e-19;
            
            wl = h * c / (this.energy * 1000 * ev2j) * 1e10;
            
        end
        
        function q = getQ(this)
            
            q = this.data.q + this.qoff;
            
        end
        
        function q = q(this)
            
            q = this.data.q + this.qoff;
            
        end
        
    end
    
end