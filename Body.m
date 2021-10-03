classdef Body < matlab.mixin.SetGet
    properties
        name
        id
        Cabs
        CR
        com
        density
        length
        mass
        meshfile
        convexmeshfile
        position
        position_w = [0;0;0]
        position_w_profile
        pos_attachments
        pos_attachments_w
        rotation
        scale
    end
    methods
        function obj = Body(Bname)
            obj.name = Bname;
        end
    end       
end