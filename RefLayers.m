classdef RefLayers < handle
    
    properties
        
        ed = [0, 0.335]
        thickness = [Inf, Inf]
        energy
        
        qoff = 0
        
    end
    
    methods
        
        function this = RefLayers(energy)
            
            this.energy = energy;
            
        end
        
        function qc = qc(this)
            
            ro = 2.818*10^-5; % radius of electron in angstrom
            k = 2 * pi / this.wavelength();
            delta = 2 * pi * ro * this.ed(end) / k^2;
            qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        function wl = wavelength(this)
            
            c = 299792458;
            h = 6.62607004e-34;
            ev2j = 1.60218e-19;
            
            wl = h * c / (this.energy * 1000 * ev2j) * 1e10;
            
        end
        
        function r = getFresnel(this, q)
            
            qz = q - this.qoff;
            r_fres = (qz - sqrt(qz.^2 - ones(size(qz)) * this.qc^2)) ./ (qz + sqrt(qz.^2 - ones(size(qz)) * this.qc^2));
            r = double(r_fres .* conj(r_fres));
            
        end
        
    end
    
end