classdef XeRefData < handle
    
    properties
        
        data
        layers
        
        qoff = 0
        qcut = 0.026
        
    end
    
    methods
        
        function this = XeRefData(file)
            
            this.data = RefData(file);
            this.layers = RefLayers(this.data);
            
        end
        
    end
    
    methods(Static)
        
        
        
    end
    
end