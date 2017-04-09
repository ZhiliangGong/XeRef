classdef RefRawData < handle
    
    properties
        
        q
        ref
        err
        energy = 10
        path
        file
        qcut = 0.026
        
    end
    
    methods
        
        function this = RefRawData(file)
            
            rawdata = importdata(file);
            this.q = rawdata(:, 1);
            this.ref = rawdata(:, 2);
            this.err = rawdata(:, 3);
            [pathname, filename, extension] = fileparts(file);
            this.file = [filename, extension];
            this.path = pathname;
            
        end
        
        function d = goodData(this)
            
            sel = this.q > this.qcut;
            d.q = this.q(sel);
            d.ref = this.ref(sel);
            d.err = this.err(sel);
            
        end
        
    end
    
end