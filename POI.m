classdef POI < dlnode
    
    properties
        w;      % wavelength
        t;      % time
        Color;  % color (dah!!)
    end
    
    methods
        
        function self = POI(w,t,color)
            self.w = w;
            self.t = t;
            self.Color = color;
        end
        
    end

end

