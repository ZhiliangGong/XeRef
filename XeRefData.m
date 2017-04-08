classdef XeRefData < handle
    
    properties
        
        rawdata
        layers
        
        qoff = 0
        qcut = 0.026
        data
        
    end
    
    methods
        
        function this = XeRefData(file)
            
            this.rawdata = RefRawData(file);
            this.layers = RefLayers(this.rawdata.energy);
            
        end
        
        function processRawData(this)
            
            sel = this.rawdata.q > this.qcut;
            this.data.q = this.rawdata.q(sel) - this.qoff;
            
        end
        
    end
    
    methods(Static)
        
        function r = getFresnelReflectivity(q, qc, qoff)
            
            qz = q - qoff;
            r_fres = (qz - sqrt(qz.^2 - ones(size(qz)) * qc^2)) ./ (qz + sqrt(qz.^2 - ones(size(qz)) * qc^2));
            r = double(r_fres .* conj(r_fres));
            
        end
        
    end
    
end