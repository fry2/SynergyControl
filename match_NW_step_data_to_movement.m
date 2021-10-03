function [outData,objNum,outDataRawAct] = match_NW_step_data_to_movement(obj,inData)
    % For inData (cell) of mean muscle activation over stride, match the swing and stance motions to the joint motion of the object.
    % Works with C:\Users\fry16\OneDrive\Documents\NW_data_local\nw_step_data.mat
    % Can be easily modified to work with normStepCell from [stepCell,normStepCell] = perStepActivation(rat,to_plot,muscNum) if nw_step_data is missing
    inRows = inData.row_names;
    inData = inData.data;
    
    % Step indices aligning with toe off and touch down. Must include end indices.
    translocs = [1,654,1034,1676,2057,2696,3064];
    
    % What are the full names for the muscle acronyms in inData? (in order)
    muscle_names = {'Vastus Lateralis';'Tibialis Anterior';'Medial Gastrocnemius';...
    'Rectus Femoris';'Biceps Femoris Posterior';'Vastus Medialis';'Semimembranosus';...
    'Illiopsoas';'Gemellus Superior';'Vastus Intermedius';'Gracilis Posterior';'Semitendinosus'};

    % Which row in inData contains the mean activation data?
    meanRow = 4;
    
    % Align the activation profiles to the object step pattern
    outData = zeros(length(obj.theta_motion),size(inData,2));
    for nwmnum = 1:size(inData,2)
        objNum(nwmnum) = obj.find_muscle_index(strrep(lower(muscle_names{nwmnum}),' ',''));
        swing = inData{meanRow,nwmnum}(1:26); stance = inData{meanRow,nwmnum}(27:end);
        full_walk_activation = [];
        for ii = 3:2:length(translocs)
            st = interp1(1:length(stance),stance,linspace(1,length(stance),translocs(ii-1)-translocs(ii-2)))';
            sw = interp1(1:length(swing),swing,linspace(1,length(swing),translocs(ii)-translocs(ii-1)))';
            full_walk_activation = [full_walk_activation;st;sw];
        end
        temp = length(obj.theta_motion)-length(full_walk_activation);
        full_walk_activation = [full_walk_activation;st(1:temp,:)];
        outData(:,nwmnum) = obj.musc_obj{objNum(nwmnum)}.ST_max.*full_walk_activation;
        outDataRawAct(:,nwmnum) = full_walk_activation;
    end
end