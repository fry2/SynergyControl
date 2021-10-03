function [obj,joint_profile,passive_tension] = design_synergy(sim_path)

    if nargin < 1
        error('Please include a simulation file path.')
    end

    if isstring(sim_path)
        sim_path = char(sim_path);
        if ~contains(sim_path,'.asim')
            % Does the simulation file path route to an ASIM file?
            error('Simulation file path invalid.')
        else
            if ~isfile(sim_path)
                % Does the routed simulation file exist?
                error('Sim file doesn''t exist.')
            end
        end
    end

    % Define Animatlab .asim file to run
    % sim_path = [fileparts(mfilename('fullpath')),'\Animatlab\','IndividualMuscleStim20190429_mod_Standalone.asim'];
    % Run .asim file with system
    sdata = processSimData(sim_path);
    aForms = {'JointMotion';'PassiveTension'};
    formInds = [find(contains({sdata.name},aForms{1})),find(contains({sdata.name},aForms{2}))];
%     if length(sdata(formInds(1)).time)~=length(sdata(formInds(2)).time)
%         error('The end time for Joint Motion and Passive Tension are different. Make sure they have the same end time.')
%     end
    for ii = 1:length(aForms)
        tempInd = formInds(ii);
        if isempty(sdata(tempInd).data)
            error(['The ',aForms{ii},' data didn''t come through in the simulation.'])
        else
            data = sdata(tempInd).data;
        end
        col_head_slice = sdata(tempInd).data_cols;
        switch aForms{ii}
            case 'JointMotion'
                    hipInd = find(contains(col_head_slice,'Hip'),1,'first');
                    kneeInd = find(contains(col_head_slice,'Knee'),1,'first');
                    ankleInd = find(contains(col_head_slice,'Ankle'),1,'first');
                joint_profile = sdata(tempInd).data;
                    % Determine which columns correspond to which joint. Animatlab doesn't necessarily follow a hip->knee->ankle convention
                    hip = joint_profile(:,hipInd);
                    knee = joint_profile(:,kneeInd);
                    ankle = joint_profile(:,ankleInd);
                joint_profile = [sdata(tempInd).time,hip,knee,ankle];
            case 'PassiveTension'
                passive_tension = cell(length(col_head_slice),2);
                for jj = 1:length(col_head_slice)
                    passive_tension{jj,1} = col_head_slice{jj};
                    passive_tension{jj,2} = data(:,jj);
                end
        end
    end
    obj = FullLeg(sim_path,joint_profile,passive_tension);
end