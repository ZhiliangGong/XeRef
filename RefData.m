classdef RefData < handle
    
    properties
        
        raw
        q
        ref
        err
        energy = 10
        qcut = 0.03
        bufferEd = 0.335
        path
        file
        
    end
    
    methods
        
        function this = RefData(file)
            
            rawdata = importdata(file);
            
            this.raw.q = rawdata(:, 1)';
            this.raw.ref = rawdata(:, 2)';
            this.raw.err = rawdata(:, 3)';
            [pathname, filename, extension] = fileparts(file);
            this.file = [filename, extension];
            this.path = pathname;
            
            sel = this.raw.q > this.qcut;
            this.q = this.raw.q(sel);
            this.ref = this.raw.ref(sel);
            this.err = this.raw.err(sel);
            
        end
        
        function d = goodData(this)
            
            sel = this.q > this.qcut;
            d.q = this.q(sel);
            d.ref = this.ref(sel);
            d.err = this.err(sel);
            
        end
        
        function d = getFND(this)
            
            fresnel = this.getFresnel();
            d.ref = this.ref ./ fresnel;
            d.err = this.err ./ fresnel;
            
        end
        
        function R = getFresnel(this)
            
            qc = this.getQc();
            
            r = (this.q - sqrt(this.q.^2 - ones(size(this.q)) * qc^2)) ./ (this.q + sqrt(this.q.^2 - ones(size(this.q)) * qc^2));
            R = r .* conj(r);
            
        end
        
        function qc = getQc(this)
            
            ro = 2.818*10^-5;
            k = 2 * pi / this.getWavelength();
            delta = 2 * pi * ro * this.bufferEd / k^2;
            qc = 2 * k * sqrt( 2 * delta );
            
        end
        
        function wl = getWavelength(this)
            
            c = 299792458;
            h = 6.62607004e-34;
            ev2j = 1.60218e-19;
            
            wl = h * c / (this.energy * 1000 * ev2j) * 1e10;
            
        end
        
    end
    
end