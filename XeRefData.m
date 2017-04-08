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
            this.layers = RefLayers(this.rawdata);
            
        end
        
    end
    
    methods(Static)
        
        
        
    end
    
end