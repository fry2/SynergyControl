classdef Muscle < matlab.mixin.SetGet
    properties
      muscle_name
      muscle_index
      enabled
      pos_attachments
      pos_attachments_w
      x_off %X offset in V (-.04=-40mV)
      max_force % Stimulus-Tension Amplitude in Newtons
      ST_max %Magnitude of the ST curve, may or not be the same as Fmax
      steepness % Steepness
      y_off % Y offset in Newtons
      RestingLength % Resting length in meters
      Kse % Serial stiffness in N/m
      Kpe % Parallel stiffness in N/m
      damping % Damping constant in Ns/m
      l_min %Minimum muscle length
      l_max %Maximum muscle length
      l_width
      mass
      opt_fiber_length
      pennation_angle
      Po
      tendon_sl
      opt_muscle_length
      vmax_fiber
      lf_lm
      muscle_length_profile
      muscle_velocity_profile
      passive_tension
    end
    methods
        function obj = Muscle()
            
        end
        function musc_vec = muscle_unit_vector(obj)
            for i = 2:size(obj.pos_attachments,1)
                temp_vec = obj.pos_attachments{i-1,1} - obj.pos_attachments{i,1};
                musc_vec = temp_vec/norm(temp_vec);
            end
        end
        
        function musc_len = muscle_length(obj)
            musc_len = 0;
            for i = 2:size(obj.pos_attachments,1)
                musc_len = musc_len + norm(obj.pos_attachments{i-1,1} - obj.pos_attachments{i,1});
            end
        end 
    end
end