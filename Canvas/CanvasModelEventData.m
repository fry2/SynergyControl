classdef (ConstructOnLoad) CanvasModelEventData < event.EventData
    
   properties
      Index
      Type
      Ends = [];
   end
   
   methods
      function data = CanvasModelEventData(type, ind, ends)
        data.Index = ind;
        data.Type = type;
        if nargin == 3
            data.Ends = ends;
        end
      end
   end
end