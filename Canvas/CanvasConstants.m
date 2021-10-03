classdef CanvasConstants < handle
        
    properties (Constant)
        %Canvas Properties
         CANVAS_LIMITS = [738 620];
        DATA_RELATIVE_PATH = fullfile('..', '..', 'data')
%         BACKGROUND_IMAGE = 'BackGround.jpg'
        BACKGROUND_IMAGE = 'White_Background.png'
        CURRENT_CANVAS_SETUP = 'Sites'
        
        %Neuron Properties
        NEURON_size = [50 50]
        NEURON_enabled = 'True';
        NEURON_restingpotential = -60;
        NEURON_timeconstant = 5;
        NEURON_initialthreshold = 50;
        NEURON_relativeaccomodation = 0.3;
        NEURON_accomodationtimeconstant = 10;
        NEURON_AHPconductance = 0;
        NEURON_AHPtimeconstant = 3;
        NEURON_ca_act_ID
        NEURON_ca_act_midpoint = -30;
        NEURON_ca_act_slope = 0.1;
        NEURON_ca_act_timeconstant = 20;
        NEURON_ca_deact_ID
        NEURON_ca_deact_midpoint = -90;
        NEURON_ca_deact_slope = -0.1;
        NEURON_ca_deact_timeconstant = 7500;
        NEURON_tonicstimulus = 0;
        NEURON_tonicnoise = 0;
        
        MUSCLE_size = [25 116.25];
        MUSCLE_transp = 'True';
        
        ADAPTER_size = [26 26];
        ADAPTER_gain_profile_ID = 'gainID';
        ADAPTER_type = 'adapter';
        
        OBSTACLE_SIZE = [25 25]
    end
    
    methods (Static)
        
%         function pos = ConstrainedPosition(type, pos)
%             if type == 'n' %target
%                 sz = CanvasConstants.NEURON_size;
%             else %obstacle
%                 sz = CanvasConstants.OBSTACLE_SIZE;
%             end
%             
%             bottomLeftLimits = [1 1] + sz/2;
%             topRightLimits = CanvasConstants.CANVAS_LIMITS - [1 1] - sz/2;
%             
%             blCond = pos < bottomLeftLimits;
%             pos(blCond) = bottomLeftLimits(blCond);
%             
%             urCond = pos > topRightLimits;
%             pos(urCond) = topRightLimits(urCond);
%             pos = round(pos);
%         end
       
    end
    
end