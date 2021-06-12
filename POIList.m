classdef POIList < dllist
 
    properties
    end
    
    methods
        
        function self = POIList()
        end
        
        function addPOI(self, w, t)
            c = self.setColor(self.Num_nodes);
            p = POI(w,t,c);
            self.add_node(p);
        end
        
        function c = setColor(self, numInstances)
            numColors = 15;
            colors = flipud(jet(numColors));
            nextColorInd = mod(numInstances, numColors)+1;
            c = colors(nextColorInd,:);
        end
        
    end
end

