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
            this.file = [filename, '.', extension];
            this.path = pathname;
            
        end
        
        function r = getFresnelReflectivity(q, qc, qoff)
            
            qz = q - qoff;
            r_fres = (qz - sqrt(qz.^2 - ones(size(qz)) * qc^2)) ./ (qz + sqrt(qz.^2 - ones(size(qz)) * qc^2));
            r = double(r_fres .* conj(r_fres));
            
        end
        
    end
    
end