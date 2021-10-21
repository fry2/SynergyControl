classdef FullLeg < matlab.mixin.SetGet
   properties
        proj_file
        original_text
        is_sim_file
        organism_name = 'Organism_1';
        organism_position;
        body_obj
        joint_obj; %cell array of joint objects (for each joint in the leg)
        musc_obj = [];
    
        torque_motion;
        theta_motion;
        theta_motion_time;
        theta_dot_motion;
        dt_motion;
        
        sampling_vector;
        passive_tension;
   end
   methods
        function obj = FullLeg(inFilePath,joint_profile,passive_tension)
        %% Process the Input File
            % Input: inFilePath: char: path to file (.aproj or .asim)
            if contains(inFilePath,'.asim')
                is_sim = 1;
                obj.is_sim_file = 1;
            elseif contains(inFilePath,'.aproj')
                is_sim = 0;
                obj.is_sim_file = 0;
            else
                error('Input inFilePath in FullLeg is not .asim or .aproj.')
            end
            
            try 
                obj.original_text = importdata(inFilePath);
                ot = obj.original_text;
            catch
                disp('No simulation file exists. Open AnimatLab and export a standalone simulation.')
                keyboard
            end
            jointTemp = ot(find(contains(ot,'<Joint>'))+1);
            joints = cell(1,length(jointTemp)+1);
            for ii = 1:length(jointTemp)
                joints{ii+1} = char(extractBetween(string(jointTemp{ii}),'>','<'));
            end
            
            if is_sim
                bodyTemp = ot(find(contains(ot,'<Type>Mesh</Type>'))-17);
            else
                bodyTemp = ot(find(contains(ot,'<Type>Mesh</Type>'))-28);
            end
            bodies = cell(1,length(bodyTemp)/2);
            counter = 1;
            for ii = 1:length(bodyTemp)
                if ~contains(lower(bodyTemp{ii}),'graphics')
                    bodies{counter} = char(extractBetween(string(bodyTemp{ii}),'>','<'));
                    counter = counter +1;
                end
            end
            %length(obj.body_obj) = length(bodies)           
        %% Initialize Body and Joint Positions in the Local Reference Frames
            %Cumulative body rotation
            cum_body_rot = eye(3);
            %Position of the joints and bodies, in their local frames.
            %obj.pos_joints = zeros(3,length(obj.body_obj));
            euler_angs_bodies = zeros(3,length(obj.body_obj));
            
            %Initialize other object values
            for ii = 1:length(bodies)
                obj.body_obj{ii,1} = Body(bodies{ii});
            end
            
            joints = joints(~cellfun(@isempty,joints));
            for ii = 1:length(joints)
                obj.joint_obj{ii,1} = JointSyn(joints{ii});
            end         
        %% Fill in Body and Joint Details
            orgIDInd = find(contains(ot,'<Organism>'),1,'first');
            if is_sim
                obj.organism_position = cellfun(@str2double,extractBetween(ot{orgIDInd+15},'"','"'))';
            elseif ~is_sim
                lpPos = find(contains(ot(orgIDInd:end),'LocalPosition'),1,'first')+orgIDInd-1;
                for tt = 1:3
                    obj.organism_position(tt) = str2double(extractBetween(string(ot{lpPos+tt}),'Actual="','"'));
                end
            end

            for i = 1:length(obj.body_obj)
                %Find the body of interest
                if ~isempty(obj.body_obj{i}.name)
                    next_body_ind = find(contains(ot,['<Name>',obj.body_obj{i}.name,'</Name>']),1,'first');
                    chain_lower_limit = next_body_ind;

                    if isempty(chain_lower_limit)
                        disp('An improper body or joint was entered, or perhaps out of order.')
                    end
                    
                    if is_sim
                        body_params = {'<ID>';'<Position ';'<Rotation ';'<Density>';'<Mass>';'<COM ';'<MeshFile>';'<ConvexMeshFile>';'<Scale '};
                    elseif ~is_sim
                        body_params = {'<ID>';'<LocalPosition>';'<Rotation>';'<Density ';'<Mass ';'<COM>';'<MeshFile>';'<ConvexMeshFile>';'<Scale>'};
                        body_param_names = {'id';'position';'rotation';'density';'mass';'com';'meshfile';'convexmeshfile';'scale'};
                    end
                    
                    for ii = 1:length(body_params)
                        pTemp = lower(erase(body_params{ii},{'<','>',' '}));
                        pInd = find(contains(ot(chain_lower_limit:end),body_params{ii}),1,'first') + chain_lower_limit - 1;
                        quoteLocs = strfind(ot{pInd},'"');
                        if ~isempty(quoteLocs)
                            sTemp = ot{pInd};
                            qTemp = reshape(quoteLocs,[2 3])';
                            parDbl = zeros(1,length(quoteLocs)/2);
                            for jj = 1:size(qTemp,1)
                                parDbl(jj) = str2double(sTemp(qTemp(jj,1)+1:qTemp(jj,2)-1));
                            end
                            obj.body_obj{i}.(pTemp) = parDbl;
                        else
                            if sum(ot{pInd}=='>')==1 && sum(ot{pInd}=='<')==1
                                parTemp = zeros(1,3);
                                for tt = 1:3
                                    parTemp(tt) = str2double(extractBetween(string(ot{pInd+tt}),'Actual="','"'));
                                end
                                if strcmp(body_param_names{ii},'rotation')
                                    obj.body_obj{i}.(body_param_names{ii}) = parTemp.*(pi/180);
                                else
                                    obj.body_obj{i}.(body_param_names{ii}) = parTemp;
                                end
                            else
                                sTemp = extractBetween(string(ot{pInd}),'>','<');
                                if isnan(str2double(sTemp))
                                    obj.body_obj{i}.(pTemp) = char(sTemp);
                                else
                                    obj.body_obj{i}.(pTemp) = str2double(sTemp);
                                end
                            end
                        end
                    end
                    
                    if length(obj.body_obj{i}.density) > 1
                        obj.body_obj{i}.density = obj.body_obj{i}.density(1);
                    end
                    
                    if length(obj.body_obj{i}.mass) > 1
                        obj.body_obj{i}.mass = obj.body_obj{i}.mass(1);
                    end
                    
                    %Find the string that lists its rotation
                    body_rot_ind = find(contains(ot(chain_lower_limit:end),'<Rotation'),1,'first') + chain_lower_limit - 1;
                    body_rot_str = ot(body_rot_ind);

                    %Read that body's rotation angles and construct a rotation matrix.
%                     quote_locs = cell2mat(strfind(body_rot_str,'"'));
%                     for j=1:length(quote_locs)/2
%                         euler_angs_bodies(j,i) = str2double(body_rot_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
%                     end
                    euler_angs_bodies(:,i) = obj.body_obj{i}.rotation';
                    %obj.euler_angs_bodies(:,i) = euler_angs_bodies(:,i);
                    %Compute the rotation matrix for each axis.
                    rotx_b = [1,0,0;0,cos(euler_angs_bodies(1,i)),-sin(euler_angs_bodies(1,i));0,sin(euler_angs_bodies(1,i)),cos(euler_angs_bodies(1,i))];
                    roty_b = [cos(euler_angs_bodies(2,i)),0,sin(euler_angs_bodies(2,i));0,1,0;-sin(euler_angs_bodies(2,i)),0,cos(euler_angs_bodies(2,i))];
                    rotz_b = [cos(euler_angs_bodies(3,i)),-sin(euler_angs_bodies(3,i)),0;sin(euler_angs_bodies(3,i)),cos(euler_angs_bodies(3,i)),0;0,0,1];

                    %Record the relative rotation of this frame with regard to the
                    %previous, as well as the absolute rotation with respect to the ground.
                    %obj.CR_bodies(:,:,i) = rotx_b*roty_b*rotz_b;
                    obj.body_obj{i}.CR = rotx_b*roty_b*rotz_b;
                    %cum_body_rot = cum_body_rot * obj.CR_bodies(:,:,i);
                    cum_body_rot = cum_body_rot * obj.body_obj{i}.CR;
%                     obj.body_obj{i) = cum_body_rot;
                    obj.body_obj{i}.Cabs = cum_body_rot;
                   
                    % Now for the joint
                    if i <= length(obj.joint_obj)
                        %Move up the minimum line value now that that body is done.
                        joint_found = contains(ot,['<Name>',obj.joint_obj{i}.name,'</Name>']);
                        %next_joint_ind = find(~cellfun(@isempty,joint_found),1,'first');
                        %chain_lower_limit = next_joint_ind;
                        obj.joint_obj{i}.index = find(joint_found,1,'first');
                        chain_lower_limit = obj.joint_obj{i}.index;

                        %Find the orientation of that joint
                        joint_rot_found = contains(ot(chain_lower_limit:end),'<Rotation');
                        joint_rot_ind = find(joint_rot_found,1,'first') + chain_lower_limit - 1;
                        joint_rot_str = ot(joint_rot_ind);
                        
                        if is_sim
                            quote_locs = cell2mat(strfind(joint_rot_str,'"'));
                            for j=1:length(quote_locs)/2
                                %obj.euler_angs_joints(j,i) = str2double(joint_rot_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                                obj.joint_obj{i}.euler_angs(j,1) = str2double(joint_rot_str{1}(quote_locs(2*j-1)+1:quote_locs(2*j)-1));
                            end
                        elseif ~is_sim
                            for tt = 1:3
                                obj.joint_obj{i}.euler_angs(tt) = (pi/180)*str2double(extractBetween(string(ot{joint_rot_ind+tt}),'Actual="','"'));
                            end
                        end
                            
                        %Find the limits of that joint             
                        if is_sim
                            obj.joint_obj{i}.limits(1,1) = double(extractBetween(string(ot{find(contains(ot(chain_lower_limit:end),'<LowerLimit>'),1,'first')+chain_lower_limit+2}),'>','<'));
                            obj.joint_obj{i}.limits(2,1) = double(extractBetween(string(ot{find(contains(ot(chain_lower_limit:end),'<UpperLimit>'),1,'first')+chain_lower_limit+2}),'>','<'));
                        elseif ~is_sim
                            obj.joint_obj{i}.limits(1,1) = (pi/180)*double(extractBetween(string(ot{find(contains(ot(chain_lower_limit:end),'<LowerLimit>'),1,'first')+chain_lower_limit+2}),'Actual="','"'));
                            obj.joint_obj{i}.limits(2,1) = (pi/180)*double(extractBetween(string(ot{find(contains(ot(chain_lower_limit:end),'<UpperLimit>'),1,'first')+chain_lower_limit+2}),'Actual="','"'));
                        end 
                        
                        % Determine if joint limits are enabled
                        joint_lim_enable = find(contains(ot(chain_lower_limit:end),'<EnableLimits'),1,'first')+chain_lower_limit-1;
                        if contains(ot{joint_lim_enable},'True')
                            obj.joint_obj{i}.enable_limit = 1;
                        end
                        
                        %Compute the rotation of the joint's axis within the local frame.
                        rotx_j = [1,0,0;0,cos(obj.joint_obj{i}.euler_angs(1)),-sin(obj.joint_obj{i}.euler_angs(1));0,sin(obj.joint_obj{i}.euler_angs(1)),cos(obj.joint_obj{i}.euler_angs(1))];
                        roty_j = [cos(obj.joint_obj{i}.euler_angs(2)),0,sin(obj.joint_obj{i}.euler_angs(2));0,1,0;-sin(obj.joint_obj{i}.euler_angs(2)),0,cos(obj.joint_obj{i}.euler_angs(2))];
                        rotz_j = [cos(obj.joint_obj{i}.euler_angs(3)),-sin(obj.joint_obj{i}.euler_angs(3)),0;sin(obj.joint_obj{i}.euler_angs(3)),cos(obj.joint_obj{i}.euler_angs(3)),0;0,0,1];

                        %Quick side track: find what kind of joint this is. If it is
                        %prismatic, it does not rotate anything. If it is a hinge, then we
                        %need to account for that rotation.
                        joint_type_found = contains(ot(chain_lower_limit:end),'<PartType');
                        joint_type_ind = find(joint_type_found,1,'first') + chain_lower_limit - 1;
                        joint_type_str = ot(joint_type_ind);

                        %pull out the joint type
                        type_begin = strfind(joint_type_str,'.');
                        type_end = strfind(joint_type_str,'<');
                        obj.joint_obj{i}.type = joint_type_str{:}(type_begin{:}(end)+1:type_end{:}(end)-1);

                        %Calculate the orientation of the axis in space, Cabs_joints, and
                        %the rotation of the joint, C_joints
                        obj.joint_obj{i}.Cabs = obj.body_obj{i}.Cabs*rotx_j*roty_j*rotz_j;
                        obj.joint_obj{i}.CR = rotx_j*roty_j*rotz_j;

                        %Find the location of the joint within the frame
                        if is_sim
                            joint_pos_str = ot(find(contains(ot(chain_lower_limit:end),'<Position'),1,'first') + chain_lower_limit - 1);
                            obj.joint_obj{i}.init_pos = cellfun(@str2double,extractBetween(joint_pos_str,'"','"'));
                        else
                            lpInd = find(contains(ot(chain_lower_limit:end),'LocalPosition'),1,'first')+chain_lower_limit-1;
                            for tt = 1:3
                                obj.joint_obj{i}.init_pos(tt,1) = str2double(extractBetween(string(ot{lpInd+tt}),'Actual="','"'));
                            end
                        end
                    end
                end
            end
            
            %%Store the world position for each body and joint
            obj.body_obj{1}.position_w = obj.organism_position';

            for i=2:length(obj.body_obj)
                %pos_bodies should be counted from 2 to end. Each column is that body's
                %position relative to the first body (which is 0s for the first).
                %Second body's position is given in the frame of the first.
                obj.body_obj{i}.position_w = obj.body_obj{i-1}.position_w + obj.body_obj{i-1}.Cabs*obj.body_obj{i}.position';
                obj.joint_obj{i-1}.init_pos_w = obj.body_obj{i}.position_w + obj.body_obj{i}.Cabs*obj.joint_obj{i-1}.init_pos;
            end
            
            % find the toe position. We approxmiate this as the distal end of the EDL. first, find the index of the EDL attachment
            edlInd = find(contains(ot,'<Name>ExtensorDigitorumLongusIns Foot LH_AnkleZ Flx</Name>'));
            % find the index of the position line for the attachment
            if obj.is_sim_file
                toePosInd = find(contains(ot(edlInd:end),'<Position'),1,'first')+edlInd-1;
                % extract the position, convert from string to double
                toePosLocal = cellfun(@str2double,extractBetween(ot{toePosInd},'"','"'));
            else
                toePosInd = find(contains(ot(edlInd:end),'<LocalPosition'),1,'first')+edlInd-1;
                for ii = 1:3
                    toePosLocal(ii,1) = double(extractBetween(string(ot{toePosInd+ii}),'Actual="','"/>'));
                end
            end
            
            initToeW = att_pos_on_demand(obj,[0,0,0],[{toePosLocal},{[]},{4},{[]}]);
            for i = 2:length(obj.body_obj)
                if i == length(obj.body_obj)
                    obj.body_obj{i}.length = norm(initToeW-obj.joint_obj{i-1}.init_pos_w);
                else
                    obj.body_obj{i}.length = norm(obj.joint_obj{i}.init_pos_w-obj.joint_obj{i-1}.init_pos_w);
                end
            end
        %% Find Muscles, Attachments, and Their Associated Bodies/Joints
            musc_found = strfind(obj.original_text,'<Type>LinearHillMuscle</Type>');
            musc_inds = find(~cellfun(@isempty,musc_found));
            musc_name_inds = musc_inds-2;

            for i = 1:length(musc_inds)
                obj.musc_obj{i,1} = Muscle();
            end
            
            name_true = 1;
            for i=1:length(obj.joint_obj)-1
                if ~strcmp(obj.joint_obj{i}.name(1:2),obj.joint_obj{i+1}.name(1:2))
                    name_true = 0;
                end
            end
            
            if name_true
                leg_name = obj.joint_obj{1}.name(1:2);
            else
                leg_name = input('What is the prefix of objects belonging to this leg?\nA cell array of prefixes may be entered. ');
            end
            
            if ischar(leg_name)
                musc_for_this_leg = strfind(obj.original_text(musc_name_inds),['<Name>',leg_name]);
            elseif iscell(leg_name)
                musc_for_this_leg = cell(length(musc_name_inds),1);
                for i=1:length(leg_name)
                    temp_musc_for_this_leg = strfind(obj.original_text(musc_name_inds),['<Name>',leg_name{i}]);
                    for j=1:length(temp_musc_for_this_leg)
                        if ~isempty(temp_musc_for_this_leg{j})
                            musc_for_this_leg{j} = 1;
                        end
                    end
                end
            end
            musc_for_this_leg_inds = cellfun(@isempty,musc_for_this_leg) == 0;
            
            %These are the indices of the names of muscles.
            %obj.musc_inds = musc_name_inds(musc_for_this_leg_inds);
            attach_to_find = cell(length(musc_inds),1);
            
            %Now that we know where the muscles are saved, we need to extract all of the attachment IDs, find their names (to figure
            %out which body they belong to), and save their locations to create joint objects.
            for i=1:length(musc_inds)
                attachment_start = contains(obj.original_text(musc_inds(i):end),'<Attachments>');
                attachment_start_ind = find(attachment_start,1,'first') + musc_inds(i) - 1;
                attachment_end = contains(obj.original_text(musc_inds(i):end),'</Attachments>');
                attachment_end_ind = find(attachment_end,1,'first') + musc_inds(i) - 1;

                attach_to_find{i} = obj.original_text(attachment_start_ind + 1:attachment_end_ind - 1);
                %remove the word "attach" from the strings
                for j=1:length(attach_to_find{i})
                    attach_to_find{i}{j} = strrep(attach_to_find{i}{j},'Attach','');
                end                
            end
            
            attach_names = cell(size(attach_to_find));
            attach_names_str = cell(size(attach_to_find));
            attach_locs = cell(size(attach_to_find));
            %find indices of the names of these attachments.
            for i=1:length(musc_inds)
                for j=1:length(attach_to_find{i})
                    %save the names
                    id_ind = find(contains(obj.original_text,strtrim(attach_to_find{i}{j})),1,'first');
                    attach_names{i}{j} = id_ind - 1;
                    attach_names_str{i}{j} = obj.original_text{attach_names{i}{j}};
                     
                    %Find the position of that attachment. This mapping is
                    %identical to the names.
                    %NEW
                    if is_sim
                        attach_pos_found = contains(obj.original_text(id_ind:end),'<Position');
                        attach_pos_ind = find(attach_pos_found,1,'first') + id_ind - 1;
                        attach_pos_str = obj.original_text(attach_pos_ind);
                        quote_locs = cell2mat(strfind(attach_pos_str,'"'));
                        for k=1:length(quote_locs)/2
                            attach_locs{i}{j}(k,1) = str2double(attach_pos_str{1}(quote_locs(2*k-1)+1:quote_locs(2*k)-1));
                        end
                    elseif ~is_sim
                        attach_pos_ind = find(contains(obj.original_text(id_ind:end),'<LocalPosition'),1,'first')+id_ind-1;
                        for tt = 1:3
                            attach_locs{i}{j}(tt,1) = str2double(extractBetween(string(ot{attach_pos_ind+tt}),'Actual="','"')); 
                        end
                    end

                end
            end
            %We will need to remove redundant attachments. If multiple muscles use the same
            %attachment, it will show up in our list twice. Therefore we
            %will start a list used_indices, and every time we add an
            %attachment's position to our record, we add the index of its
            %name to this list. Attachments whose indices already appear on
            %the list will be ignored.
            used_indices = [];
            
            %Tally up the total number of attachments on each body so that
            %we can set the size for pos_attachments
            %Cells don't allow you to increment the size
            pelvis_count = 0;
            femur_count = 0;
            tibia_count = 0;
            foot_count = 0;
            for j=1:length(attach_names)
                for k=1:length(attach_names{j})
                    if isempty(find(used_indices == attach_names{j}{k},1,'first'))
                        temp_name_str = char(obj.original_text(attach_names{j}{k}));
                        if contains(temp_name_str,' Pelvis ')
                            pelvis_count = pelvis_count + 1;
                        elseif contains(temp_name_str,' Femur ')
                            femur_count = femur_count + 1;
                        elseif contains(temp_name_str,' Tibia ')
                            %spaces around Tibia important bc muscle
                            %'Tibialis Anterior' contains Tibia
                            tibia_count = tibia_count + 1;
                        elseif contains(temp_name_str,' Foot ')
                            foot_count = foot_count + 1;
                        else
                            keyboard
                            error('Attachment not labeled with a body part')
                        end
                        used_indices = [used_indices,attach_names{j}{k}];
                    end
                end
            end
            used_indices = used_indices';
            saved_indices = used_indices';
            %Set the size of pos_attachments for each body
            obj.body_obj{1}.pos_attachments = cell(pelvis_count,3);
            obj.body_obj{2}.pos_attachments = cell(femur_count,3);
            obj.body_obj{3}.pos_attachments = cell(tibia_count,3);
            obj.body_obj{4}.pos_attachments = cell(foot_count,3);
            
            %Store information about attachment points for each body in
            %pos_attachments. Resultant cell array holds position
            %information for XYZ and name strings of every unique
            %attachment point on each body
            for j=1:length(attach_names)
                obj.musc_obj{j}.muscle_index = musc_name_inds(j);
                for k=1:length(attach_names{j})
                    if isequal(size(attach_locs{j}{k}),[3,1])
                        temp_name_str = char(extractBetween(string(obj.original_text{attach_names{j}{k}}),'<Name>','</Name>'));
                        if contains(temp_name_str,' Pelvis ')
                            body = 2;
                        elseif contains(temp_name_str,' Femur ')
                            body = 3;
                        elseif contains(temp_name_str,' Tibia ')
                            body = 4;
                        elseif contains(temp_name_str,' Foot ')
                            body = 5;
                        else
                            keyboard
                            error('Attachment not labeled with a body part')
                        end
                        %Following if statement necessary to save
                        %information only for non-redundant attachment
                        %points
                        obj.musc_obj{j,1}.pos_attachments{k,1} = attach_locs{j}{k};
                        num_attach = strcat(num2str(size(obj.musc_obj{j,1}.pos_attachments{k,1},2)),') ');
                        obj.musc_obj{j,1}.pos_attachments{k,2} = [num_attach,obj.original_text{attach_names{j}{k}}(7:end-7),'\n'];
                        obj.musc_obj{j,1}.pos_attachments{k,3} = body-1;
                        obj.musc_obj{j,1}.pos_attachments{k,4} = [];
                        if ~isempty(find(used_indices == attach_names{j}{k},1,'first'))
                            used_indices(find(used_indices == attach_names{j}{k},1,'first')) = 0;
                            jj = sum(~cellfun(@isempty,obj.body_obj{body-1}.pos_attachments(:,1)))+1;
                            obj.body_obj{body-1}.pos_attachments{jj,1} = [obj.body_obj{body-1}.pos_attachments{jj,1},attach_locs{j}{k}];
                            num_attach = strcat(num2str(size(obj.body_obj{body-1}.pos_attachments{jj,1},2)),') ');
                            obj.body_obj{body-1}.pos_attachments{jj,2} = [obj.body_obj{body-1}.pos_attachments{jj,2},num_attach,obj.original_text{attach_names{j}{k}}(7:end-7),'\n'];
                        else
                            %skip this reused attachment
                        end
                    end  
                end
            end

            for i=1:length(obj.body_obj)
                obj.body_obj{i}.pos_attachments_w = obj.body_obj{i}.pos_attachments;
                for j = 1:size(obj.body_obj{i}.pos_attachments,1)
                    obj.body_obj{i}.pos_attachments_w{j,1} = obj.body_obj{i}.position_w + obj.body_obj{i}.Cabs*obj.body_obj{i}.pos_attachments_w{j,1};
                end
            end       
        %% Store joint properties
            tstart = tic;
            %%Neutral position. Enable if you want to have a contant profile where all joints are at 90 degrees
            if isempty(joint_profile)
                if is_sim
                    simEndTime = double(extractBetween(string(ot{contains(ot,'<SimEndTime>')}),'>','<'));
                    physTS = double(extractBetween(string(ot{contains(ot,'<PhysicsTimeStep>')}),'>','<'));
                elseif ~is_sim
                    simEndTime = double(extractBetween(string(ot{contains(ot,'<SimEndTime ')}),'Actual="','"'));
                    physTS = double(extractBetween(string(ot{contains(ot,'<PhysicsTimeStep ')}),'Actual="','"'));
                end
                obj.dt_motion = physTS;
                timeVec = 0:physTS:simEndTime-.1;
                % If no joint profile is provided, assume that the leg should be in a neutral position (halfway between limits) and stationary
                jangTemp = zeros(1,3);
                for ii = 1:3
                    jangTemp(ii) = sum(obj.joint_obj{ii}.limits)/2;
                end
                %jangTemp = [0.2793    0.4363    0.7113];
                %jangTemp = [-0.1745   -0.4363   -0.3491];
                joint_profile = [timeVec',repmat(jangTemp,length(timeVec),1)];
            else
                obj.dt_motion = joint_profile(10,1)-joint_profile(9,1);
            end

            [jointMotion,jointMotionDot] = store_joint_rotmat_profile(obj,joint_profile);
            store_torque_and_kinematic_data(obj,jointMotion,jointMotionDot);
            store_sampling_vector(obj);
            store_jointbodyw_position_profiles(obj,jointMotion);
            store_joint_params(obj);
            telapsed = toc(tstart);
            %disp(['Joint Properties Stored.',' (',num2str(telapsed),'s)'])
            clear tstart telapsed
        %% Store muscle properties  
            tstart = tic;
            store_animatlab_params(obj);
            store_johnson_params(obj);
            store_muscle_profiles(obj);
            store_input_muscle_passive_tension(obj,passive_tension);
            telapsed = toc(tstart);
            %disp(['Muscle Properties Stored.',' (',num2str(telapsed),'s)'])
            clear tstart telapsed
        end
        %% Function: Store the joint rotation matrices in the joint objects
        function [jointMotion,jointMotionDot] = store_joint_rotmat_profile(obj,jointMotion)     
            ot = obj.original_text;
            t = jointMotion(:,1);
            jointMotionDot = [jointMotion(1:length(t)-1,1),diff(jointMotion(1:length(t),2:4))/obj.dt_motion];
            jointMotion = jointMotion(1:length(t),:);
            %% For Moving a single joint at a time
            onejointonly = 0;
            if onejointonly
                load('joint_angles_maxmin.mat')
                fullrangewalking = 1;
                standardwalking = 0;
                hiponly = 0;
                kneeonly = 0;
                ankleonly = 1;
                if hiponly+kneeonly+ankleonly > 1
                    error('Only one active at a time')
                end
                if standardwalking
                    if hiponly == 1
                        joint_profile = jointMotion(200:end,1:4);
                        joint_profile(:,3) = max(joint_profile(:,3))-min(joint_profile(:,3));
                        joint_profile(:,4) = max(joint_profile(:,4))-min(joint_profile(:,4));
                        jointMotionDot = jointMotionDot(200:end,1:4);
                        jointMotionDot(:,3) = max(jointMotionDot(:,3))-min(jointMotionDot(:,3));
                        jointMotionDot(:,4) = max(jointMotionDot(:,4))-min(jointMotionDot(:,4));
                    elseif kneeonly == 1
                        joint_profile = jointMotion(200:end,1:4);
                        joint_profile(:,2) = max(joint_profile(:,2))-min(joint_profile(:,2));
                        joint_profile(:,4) = max(joint_profile(:,4))-min(joint_profile(:,4));
                        jointMotionDot = jointMotionDot(200:end,1:4);
                        jointMotionDot(:,2) = max(jointMotionDot(:,2))-min(jointMotionDot(:,2));
                        jointMotionDot(:,4) = max(jointMotionDot(:,4))-min(jointMotionDot(:,4));
                    elseif ankleonly == 1
                        joint_profile = jointMotion(200:end,1:4);
                        joint_profile(:,2) = max(joint_profile(:,2))-min(joint_profile(:,2));
                        joint_profile(:,3) = max(joint_profile(:,3))-min(joint_profile(:,3));
                        jointMotionDot = jointMotionDot(200:end,1:4);
                        jointMotionDot(:,2) = max(jointMotionDot(:,2))-min(jointMotionDot(:,2));
                        jointMotionDot(:,3) = max(jointMotionDot(:,3))-min(jointMotionDot(:,3));
                    end
                elseif fullrangewalking
                    limits = (pi/180)*[30.00 -50.00;...
                        -63.42 36.58;...
                        29.95 -70.05];
                    ranges = limits(:,1)-limits(:,2);
                    offsets = ranges/2-limits(:,1);
                    joint_profile = zeros(size(jointMotion));
                    joint_profile(:,1) = jointMotion(:,1);
                    joint_profile(:,2) = joint_profile(:,2)+(ranges(1)/2*sin(.1*2*pi*jointMotion(:,1))-offsets(1));
                    joint_profile(:,3) = joint_profile(:,3)+(ranges(2)/2*sin(.1*2*pi*jointMotion(:,1))-offsets(2));
                    joint_profile(:,4) = joint_profile(:,4)+(ranges(3)/2*sin(.1*2*pi*jointMotion(:,1))-offsets(3));
                    jointMotionDot = jointMotionDot(200:end,1:4);
                    if hiponly == 1
                        jointMotionDot = jointMotionDot(200:end,1:4);
%                         joint_profile(:,3) = max(joint_profile(:,3))-min(joint_profile(:,3));
%                         joint_profile(:,4) = max(joint_profile(:,4))-min(joint_profile(:,4));
                        %Set the knee and ankle to the C/J zero angle as defined in Animatlab
                        joint_profile(:,3) = (pi/180)*81.578;
                        joint_profile(:,4) = (pi/180)*-20.053;
                        jointMotionDot(:,3) = max(jointMotionDot(:,3))-min(jointMotionDot(:,3));
                        jointMotionDot(:,4) = max(jointMotionDot(:,4))-min(jointMotionDot(:,4));
                    elseif kneeonly == 1
                        jointMotionDot = jointMotionDot(200:end,1:4);
%                         joint_profile(:,2) = max(joint_profile(:,2))-min(joint_profile(:,2));
%                         joint_profile(:,4) = max(joint_profile(:,4))-min(joint_profile(:,4));
                        joint_profile(:,2) = 0;
                        joint_profile(:,4) =(pi/180)*81.578;
                        jointMotionDot(:,2) = max(jointMotionDot(:,2))-min(jointMotionDot(:,2));
                        jointMotionDot(:,4) = max(jointMotionDot(:,4))-min(jointMotionDot(:,4));
                    elseif ankleonly == 1
                        jointMotionDot = jointMotionDot(200:end,1:4);
%                         joint_profile(:,2) = max(joint_profile(:,2))-min(joint_profile(:,2));
%                         joint_profile(:,3) = max(joint_profile(:,3))-min(joint_profile(:,3));
                        joint_profile(:,2) = 0;
                        joint_profile(:,3) = (pi/180)*81.578;
                        jointMotionDot(:,2) = max(jointMotionDot(:,2))-min(jointMotionDot(:,2));
                        jointMotionDot(:,3) = max(jointMotionDot(:,3))-min(jointMotionDot(:,3));
                    end
                end
            end
            %% Continuing on once joint profile has been established
            for i = 1:length(obj.body_obj)-1
                %uu_joint(:,i) = obj.CR_bodies(:,:,i+1)*obj.three_axis_rotation(obj.euler_angs_joints(:,i+1))*[-1;0;0];
                %These are defined in the distal body's coordinate system (femur for the hip, etc).
                %If you want uu_joint in world coords, you need to pre-multiply by each previous body's obj.CR_bodies(:,:,i)
                %uu_joint == Ext/Flx
                %uu_joint2 == ExR/InR
                %uu_joint3 == Add/Abd
                %axismatrix = -obj.CR_bodies(:,:,i+1)*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs);
                axismatrix = -obj.body_obj{i+1}.CR*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs);
                     obj.joint_obj{i}.uu_joint  = axismatrix(:,1);
                     %Need these?
                    obj.joint_obj{i}.uu_joint2 = axismatrix(:,2);
                    obj.joint_obj{i}.uu_joint3 = axismatrix(:,3);
                    obj.joint_obj{i}.init_rot = jointMotion(1,i+1);
            end
            %This is the same thing as the above loops. The euler angs bodies term is the same as CR_bodies. This is how Alex solved for the u_joint axes
            %u_joint = obj.three_axis_rotation(obj.euler_angs_bodies(:,2))*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs)*[-1;0;0];
            a = zeros(3,3,size(jointMotion,1)); 
            b = zeros(3,3,size(jointMotion,1)); 
            c = zeros(3,3,size(jointMotion,1)); 
            hip_uu2 = zeros(3,size(jointMotion,1));
            hip_uu3 = hip_uu2;
            knee_uu2 = hip_uu2+obj.joint_obj{2}.uu_joint2;
            knee_uu3 = hip_uu2+obj.joint_obj{2}.uu_joint3;
            ankle_uu2 = hip_uu2+obj.joint_obj{3}.uu_joint2;
            ankle_uu3 = hip_uu2+obj.joint_obj{3}.uu_joint3;
            
            for i=1:size(jointMotion,1)
                a(:,:,i) = axis_angle_rotation(obj,jointMotion(i,2),obj.joint_obj{1}.uu_joint);
                b(:,:,i) = axis_angle_rotation(obj,jointMotion(i,3),obj.joint_obj{2}.uu_joint);
                c(:,:,i) = axis_angle_rotation(obj,jointMotion(i,4),obj.joint_obj{3}.uu_joint);
                %knee  = -1*b(:,:,i)*obj.CR_bodies(:,:,3)*obj.three_axis_rotation(obj.joint_obj{2}.euler_angs);
                knee  = -1*b(:,:,i)*obj.body_obj{3}.CR*obj.three_axis_rotation(obj.joint_obj{2}.euler_angs);
                %ankle = -1*obj.CR_bodies(:,:,4)*obj.three_axis_rotation(obj.joint_obj{3}.euler_angs);
                ankle = -1*obj.body_obj{4}.CR*obj.three_axis_rotation(obj.joint_obj{3}.euler_angs);
                knee_uu2(:,i)  = knee(:,2);
                knee_uu3(:,i)  = knee(:,3);
                ankle_uu2(:,i) = ankle(:,2);
                ankle_uu3(:,i) = ankle(:,3);
            end
%             hip = -1*a(:,:,1)*obj.CR_bodies(:,:,2)*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs);
%             hip = a(:,:,1)*obj.CR_bodies(:,:,2)*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs);
            hip = a(:,:,1)*obj.body_obj{2}.CR*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs);
            hip_uu2 = hip_uu2-hip(:,2);
            hip_uu3 = hip_uu3-hip(:,3);
            
            obj.joint_obj{1}.joint_rotmat_profile = a;
            obj.joint_obj{2}.joint_rotmat_profile = b;
            obj.joint_obj{3}.joint_rotmat_profile = c;
            
                obj.joint_obj{1}.uu_joint2  = hip_uu2;
                obj.joint_obj{1}.uu_joint3  = hip_uu3;
%                 obj.joint_obj{1}.uuw_joint = zeros(3,size(jointMotion,1))+obj.CR_bodies(:,:,1)*obj.joint_obj{1}.uu_joint;
                obj.joint_obj{1}.uuw_joint = zeros(3,size(jointMotion,1))+obj.body_obj{1}.CR*obj.joint_obj{1}.uu_joint;
                obj.joint_obj{1}.uuw_joint2 = obj.body_obj{1}.CR*hip_uu2;
                obj.joint_obj{1}.uuw_joint3 = obj.body_obj{1}.CR*hip_uu3;
%             kneeworld = obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2);
              kneeworld = obj.body_obj{1}.CR*obj.body_obj{2}.CR;
                obj.joint_obj{2}.uu_joint2  = knee_uu2;
                obj.joint_obj{2}.uu_joint3  = knee_uu3;
                obj.joint_obj{2}.uuw_joint  = zeros(3,size(jointMotion,1))+kneeworld*obj.joint_obj{2}.uu_joint;
                obj.joint_obj{2}.uuw_joint2 = kneeworld*knee_uu2;
                obj.joint_obj{2}.uuw_joint3 = kneeworld*knee_uu3;
            %ankleworld = obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3);
            ankleworld = obj.body_obj{1}.CR*obj.body_obj{2}.CR*obj.body_obj{3}.CR;
                obj.joint_obj{3}.uu_joint2  = ankle_uu2;
                obj.joint_obj{3}.uu_joint3  = ankle_uu3;
                obj.joint_obj{3}.uuw_joint  = zeros(3,size(jointMotion,1))+ankleworld*obj.joint_obj{3}.uu_joint;
                obj.joint_obj{3}.uuw_joint2 = ankleworld*ankle_uu2;
                obj.joint_obj{3}.uuw_joint3 = ankleworld*ankle_uu3;
        end
        %% Function: On Demand: Joint World Axes on Demand
        function axesMat = joint_axes_from_angles(obj,theta)
            warnFlag = 0;
            for ii = 1:3
                limBool = [theta(ii) > max(obj.joint_obj{ii}.limits) theta(ii) < min(obj.joint_obj{ii}.limits)];
                if any(limBool)
                    warnFlag = 1;
                    theta(ii) = obj.joint_obj{ii}.limits(limBool);
                end
            end
            
            if warnFlag
                warning('Function: FullLeg.joint_axes_from_angles: desired theta is outside joint limits.')
            end
            
%             axesMat(:,1) = obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs)*[-1;0;0];
%             axesMat(:,2) = axis_angle_rotation(obj,theta(1),axesMat(:,1))*obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.three_axis_rotation(obj.joint_obj{2}.euler_angs)*[-1;0;0];
%             axesMat(:,3) = axis_angle_rotation(obj,theta(2),axesMat(:,2))*obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.CR_bodies(:,:,4)*obj.three_axis_rotation(obj.joint_obj{3}.euler_angs)*[-1;0;0];
            axesMat(:,1) = obj.body_obj{1}.CR*obj.body_obj{2}.CR*obj.three_axis_rotation(obj.joint_obj{1}.euler_angs)*[-1;0;0];
            axesMat(:,2) = axis_angle_rotation(obj,theta(1),axesMat(:,1))*obj.body_obj{1}.CR*obj.body_obj{2}.CR*obj.body_obj{3}.CR*obj.three_axis_rotation(obj.joint_obj{2}.euler_angs)*[-1;0;0];
            axesMat(:,3) = axis_angle_rotation(obj,theta(2),axesMat(:,2))*obj.body_obj{1}.CR*obj.body_obj{2}.CR*obj.body_obj{3}.CR*obj.body_obj{4}.CR*obj.three_axis_rotation(obj.joint_obj{3}.euler_angs)*[-1;0;0];
        end
        %% Function: Store the Torque and Kinematic Data from Torque, Theta, ThetaDot input
        function store_torque_and_kinematic_data(obj,theta,theta_dot)
            %Saves torque, angle, and velocity data from desired motions
            %obj.torque_motion = torque;
            %obj.theta_motion = theta;
            %obj.theta_dot_motion = theta_dot;
            
            %Torque data
            % Torque for only stance in units of [Nmm]
                    %load('TorqueinNmm.mat');
            % Simulated torque for stance and swing in units of [Nm]
            % Generated by ConnectData.m in the Rat/Main Sim and Ani compare files
            % TorqueAll refers to torque for all joints for both stance and swing
            % TorqueinNm is -1*measured stance torque. The swing portion of torque is generated by a Simulink model
            if ismac
                load('/Volumes/GoogleDrive/My Drive/Rat/Main Sim and Ani compare files/AllTraining20180824.mat','TorqueAll');
            else
                load([fileparts(pwd),'\SynergyControl\Data\AllTraining20180824.mat'],'TorqueAll');
            end
            torque = TorqueAll;
            for i=1:length(obj.body_obj)-1
                obj.joint_obj{i}.rec_angle_profile = theta(:,i+1);
                obj.joint_obj{i}.rec_angle_time = theta(:,1);
                obj.joint_obj{i}.rec_angleDot_profile = theta_dot(:,i+1);
                obj.joint_obj{i}.rec_torque_profile = torque(:,i);
            end
            obj.theta_motion = theta(:,2:4);
            obj.theta_motion_time = theta(:,1);
            obj.theta_dot_motion = theta_dot(:,2:4);
            obj.torque_motion = TorqueAll;
        end
        %% Function: Store Joint Position Profiles in Joint Objects
        function store_jointbodyw_position_profiles(obj,joint_profile)
            sizer    = zeros(size(joint_profile,1),3);
            kneepos  = sizer;
            anklepos = sizer;
            footpos  = sizer;
            tibiapos = sizer;
            
            hippos = sizer+obj.joint_obj{1}.init_pos_w';
            femurpos = sizer+(obj.body_obj{1}.CR*obj.body_obj{2}.position')';
            for i=1:size(joint_profile,1)
                a = axis_angle_rotation(obj,joint_profile(i,2),obj.joint_obj{1}.uu_joint);
                b = axis_angle_rotation(obj,joint_profile(i,3),obj.joint_obj{2}.uu_joint);
                c = axis_angle_rotation(obj,joint_profile(i,4),obj.joint_obj{3}.uu_joint);
                
                tibiapos(i,:) = obj.body_obj{1}.CR*obj.body_obj{2}.position'+obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*obj.body_obj{3}.position';
                kneerel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR;
                kneepos(i,:) = tibiapos(i,:)+(kneerel*obj.joint_obj{2}.init_pos)';
                footpos(i,:) = tibiapos(i,:)+(kneerel*obj.body_obj{4}.position')';
                anklerel = kneerel*c*obj.body_obj{4}.CR;
                anklepos(i,:) = footpos(i,:)+ (anklerel*obj.joint_obj{3}.init_pos)';  
                %kneepos =  kneepos(i,:)*1000;
                %jointangles = joint_profile(i,2:4)';
            end
            obj.body_obj{1}.position_w_profile = sizer;
            obj.body_obj{2}.position_w_profile = femurpos;
            obj.body_obj{3}.position_w_profile = tibiapos;
            obj.body_obj{4}.position_w_profile = footpos;
            obj.joint_obj{1}.sim_position_profile = hippos;
            obj.joint_obj{2}.sim_position_profile = kneepos;
            %obj.joint_obj{2}.sim_position_profile = kneepos+.001*ones(length(kneepos),1);
            obj.joint_obj{3}.sim_position_profile = anklepos;
        end
        %% Function: On Demand: Joint Position on Demand
        function jointMat = joint_pos_on_demand(obj,theta)
            % For an input theta vector, output world positions of the three joints
            warnFlag = 0;
            for ii = 1:3
                limBool = [theta(ii) > max(obj.joint_obj{ii}.limits) theta(ii) < min(obj.joint_obj{ii}.limits)];
                if any(limBool)
                    warnFlag = 1;
                    theta(ii) = obj.joint_obj{ii}.limits(limBool);
                end
            end
            
            if warnFlag
                warning('Function: FullLeg.joint_pos_on_demand: desired theta is outside joint limits.')
            end
            
            [r,c] = size(theta);
            if ~(r==3 && c==1) && ~(r==1 && c==3)
                if any([r,c]==3)
                    if r == 3
                        theta = theta(:,1);
                    elseif c == 3
                        theta = theta(1,:);
                    else
                        error('weird error')
                    end
                else
                    error('Function: FullLeg.joint_pos_on_demand: input theta vector is not 1x3')
                end
            end
            
            axesMat = zeros(3,3);
            for i = 1:3
%                 axesMat(:,i) = obj.CR_bodies(:,:,i+1)*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
                axesMat(:,i) = obj.body_obj{i+1}.CR*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
            end
                
            a = axis_angle_rotation(obj,theta(1),axesMat(:,1));
            b = axis_angle_rotation(obj,theta(2),axesMat(:,2));
            c = axis_angle_rotation(obj,theta(3),axesMat(:,3));

            tibiapos = obj.body_obj{1}.CR*obj.body_obj{2}.position'+obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*obj.body_obj{3}.position';
            kneerel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR;
            kneepos = tibiapos+(kneerel*obj.joint_obj{2}.init_pos);
            footpos = tibiapos+(kneerel*obj.body_obj{4}.position');
            anklerel = kneerel*c*obj.body_obj{4}.CR;
            anklepos = footpos+(anklerel*obj.joint_obj{3}.init_pos);

            jointMat(:,1) = obj.joint_obj{1}.init_pos_w;
            jointMat(:,2) = kneepos;
            jointMat(:,3) = anklepos;
        end
        %% Function: On Demand: Muscle Attachment Position on Demand
        function outPos = att_pos_on_demand(obj,theta,muscAtt)
            %Input: muscAtt: 1x4 cell array with information about muscle attachment point
                % [attachment intial position, attachment name, body number of attachment, attachment profile for input waveform]
            for ii = 1:3
                limBool = [theta(ii) > max(obj.joint_obj{ii}.limits) theta(ii) < min(obj.joint_obj{ii}.limits)];
                if any(limBool) && obj.joint_obj{ii}.enable_limit
                    warning(['Desired theta ',num2str(theta(ii)),' for ',obj.joint_obj{ii}.name,' is outside joint limits, '...
                        '[',num2str(obj.joint_obj{ii}.limits(1)),', ',num2str(obj.joint_obj{ii}.limits(2)),']'])
                    theta(ii) = obj.joint_obj{ii}.limits(limBool);
                end
            end
            
            [r,c] = size(theta);
            if ~(r==3 && c==1) && ~(r==1 && c==3)
                if any([r,c]==3)
                    if r == 3
                        theta = theta(:,1);
                    elseif c == 3
                        theta = theta(1,:);
                    else
                        error('weird error')
                    end
                else
                    error('Function: FullLeg.att_pos_on_demand: input theta vector is not 1x3')
                end
            end
            
            outPos = zeros(3,1);
            axesMat = zeros(3,3); 
            bodyNum = muscAtt{1,3};

            for i = 1:3
                axesMat(:,i) = obj.body_obj{i+1}.CR*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
            end
                
            a = axis_angle_rotation(obj,theta(1),axesMat(:,1));
            b = axis_angle_rotation(obj,theta(2),axesMat(:,2));
            c = axis_angle_rotation(obj,theta(3),axesMat(:,3));
            
            pelPos = obj.organism_position';
            
            femurpos = pelPos+(obj.body_obj{1}.CR*obj.body_obj{2}.position');
            hiprel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR;
            
            tibiapos = femurpos+hiprel*obj.body_obj{3}.position';
            kneerel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR;
            
            footpos = tibiapos+(kneerel*obj.body_obj{4}.position');
            anklerel = kneerel*c*obj.body_obj{4}.CR;
            
            switch bodyNum
                case 1
                    outPos = pelPos+obj.body_obj{1}.CR*muscAtt{1};
                case 2
                    outPos = (femurpos+hiprel*muscAtt{1});
                case 3
                    outPos = (tibiapos+kneerel*muscAtt{1});
                case 4
                    outPos = footpos+(anklerel*muscAtt{1});
            end
        end
        %% Function: On Demand: Body COM Position on Demand
        function outPos = com_pos_on_demand(obj,theta,bodyNum)
            % Provides the world position of the COM for a desired body
            % Input: muscAtt: 1x4 cell array with information about muscle attachment point
                % [attachment intial position, attachment name, body number of attachment, attachment profile for input waveform]
                warnFlag = 0;
            for ii = 1:3
                limBool = [theta(ii) > max(obj.joint_obj{ii}.limits) theta(ii) < min(obj.joint_obj{ii}.limits)];
                if any(limBool)
                    warnFlag = 1;
                    theta(ii) = obj.joint_obj{ii}.limits(limBool);
                end
            end
            
            if warnFlag
                warning('Function: FullLeg.att_pos_on_demand: desired theta is outside joint limits.')
            end
            
            [r,c] = size(theta);
            if ~(r==3 && c==1) && ~(r==1 && c==3)
                if any([r,c]==3)
                    if r == 3
                        theta = theta(:,1);
                    elseif c == 3
                        theta = theta(1,:);
                    else
                        error('weird error')
                    end
                else
                    error('Function: FullLeg.att_pos_on_demand: input theta vector is not 1x3')
                end
            end
            
            outPos = zeros(3,1);
            axesMat = zeros(3,3); 
            %bodyNum = muscAtt{1,3};

            for i = 1:3
                axesMat(:,i) = obj.body_obj{i+1}.CR*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
            end
                
            a = axis_angle_rotation(obj,theta(1),axesMat(:,1));
            b = axis_angle_rotation(obj,theta(2),axesMat(:,2));
            c = axis_angle_rotation(obj,theta(3),axesMat(:,3));
            
            femurpos = obj.body_obj{1}.CR*obj.body_obj{2}.position';
            hiprel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR;
            
            tibiapos = femurpos+hiprel*obj.body_obj{3}.position';
            kneerel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR;
            %kneepos = tibiapos+(kneerel*obj.joint_obj{2}.init_pos);
            
            footpos = tibiapos+(kneerel*obj.body_obj{4}.position');
            anklerel = kneerel*c*obj.body_obj{4}.CR;
            %anklepos = footpos+(anklerel*obj.joint_obj{3}.init_pos);
            
            switch bodyNum
                case 2
                    outPos = (femurpos+hiprel*obj.body_obj{2}.com');
                case 3
                    outPos = (tibiapos+kneerel*obj.body_obj{3}.com');
                case 4
                    outPos = footpos+(anklerel*obj.body_obj{4}.com');
            end
            
        end
        %% Function: Store Animatlab Parameters Simulation File in the Muscle Objects
        function store_animatlab_params( obj )

            num_musc = length(obj.musc_obj);
            is_sim = obj.is_sim_file;
            
            load([fileparts(pwd),'\SynergyControl\Data\AnimatlabProperties.mat']);
            
            property_names = Properties{4,1}(:,1);
            numProps = length(property_names);
            property_values = zeros(num_musc*numProps,1);
            property_inds = zeros(num_musc*numProps,3);

            %for each property
            for i=1:num_musc
                %Define its line as the "lower limit." All properties will be found
                %after it in the file.
                lower_limit = obj.musc_obj{i}.muscle_index;
                for j=1:numProps
                    if strcmp(property_names{j},'damping')
                        %We should find muscle damping, which is actually called "B,"
                        %just like the muscle activation's curve. Therefore we need to
                        %be a little more specific. We want the second B that shows up
                        %for a particular muscle.
                        if is_sim
                            prop_to_find = '<B>';
                        elseif ~is_sim
                            prop_to_find = '<B Value=';
                        end
                        which_to_find = 2;
                    else
                        %Now find the desired property, formatted like the file.
                        if is_sim
                            prop_to_find = ['<',property_names{j},'>'];
                        elseif ~is_sim
                            prop_to_find = ['<',property_names{j},' Value='];
                        end
                        which_to_find = 1;
                    end

                    %Find that string in the file, the first time after the lower limit
                    %(where the named object is found). 
                    prop_found = find(contains(obj.original_text(lower_limit:end),prop_to_find),which_to_find) + lower_limit - 1;

                    %Find the index at which this occurs, and save this for quick
                    %reference later. Remember that we were only examining
                    %original_text after the lower_limit, so we need to add that back
                    %on. -1 makes the indexing work properly.
                    %temp = find(~cellfun(@isempty,prop_found)) + lower_limit - 1;
                    property_inds(numProps*(i-1)+j,1) = prop_found(which_to_find);
                    
                    if is_sim
                        %Find the final index of the row to keep before the number begins
                        %Number we're looking for is formatted like
                        %'<A>-0.04<A>'
                        %Index of > before the number we want
                        property_inds(numProps*(i-1)+j,2) = length(prop_to_find);

                        %Find the first index of the row to keep after the number begins
                        %Index of < after number we want
                        property_inds(numProps*(i-1)+j,3) = cell2mat(strfind(obj.original_text(property_inds(numProps*(i-1)+j,1)),'</'));
                    elseif ~is_sim
                        qInds = find(obj.original_text{property_inds(numProps*(i-1)+j,1)}=='"');
                        property_inds(numProps*(i-1)+j,2:3) = qInds(end-1:end);
                    end
                end
            end

            %For each property that we're changing, read in the value from the current simulation
            for i=1:size(property_inds,1) %Each property
                property_values(i) = str2double(obj.original_text{property_inds(i,1)}(property_inds(i,2)+1:property_inds(i,3)-1));
            end

            %Make a case-insensitive list of all the muscles on this leg
            
            for ii = 1:length(obj.musc_obj)
                musc_inds(ii) = obj.musc_obj{ii}.muscle_index;
            end
            
            musc_names = obj.original_text(musc_inds);
            for i=1:num_musc
                obj.musc_obj{i,1}.muscle_name = lower(musc_names{i}(7:end-7));
                obj.musc_obj{i,1}.muscle_index = musc_inds(i);
                obj.musc_obj{i,1}.enabled = contains(obj.original_text{find(contains(obj.original_text(musc_inds(i):end),'Enabled'),1,'first')+musc_inds(i)-1},'True');
                %Store the muscle properties in individual muscle objects
                obj.musc_obj{i}.x_off = property_values(numProps*i-9);
                obj.musc_obj{i}.ST_max = property_values(numProps*i-8);
                obj.musc_obj{i}.steepness = property_values(numProps*i-7);
                obj.musc_obj{i}.y_off = property_values(numProps*i-6);
                obj.musc_obj{i}.RestingLength = property_values(numProps*i-5);
                obj.musc_obj{i}.l_width = property_values(numProps*i-4);
                obj.musc_obj{i}.Kse = property_values(numProps*i-3);
                obj.musc_obj{i}.Kpe = property_values(numProps*i-2);
                obj.musc_obj{i}.damping = property_values(numProps*i-1);
                obj.musc_obj{i}.max_force = property_values(numProps*i);
            end
        end
        %% Function: Store Joint Parameters from Animatlab Sim File
        function store_joint_params(obj)
            numJoints = length(obj.joint_obj);
            ot = obj.original_text;
            for ii = 1:numJoints
                joint = obj.joint_obj{ii};
                jInd = obj.joint_obj{ii}.index;
                fricCoeff = find(contains(ot(jInd:end),'<Coefficient'),1,'first')+jInd-1;
                if obj.is_sim_file
                    joint.fricCoeff = double(extractBetween(string(ot{fricCoeff}),'>','</'));
                elseif ~obj.is_sim_file
                    joint.fricCoeff = double(extractBetween(string(ot{fricCoeff}),'Actual="','"'));
                end
                obj.joint_obj{ii} = joint;
            end
        end
        %% Function: Store Muscle Parameters from Johnson 2011 Paper in the Muscle Objects
        function store_johnson_params(obj)
            params = [];
            num_musc = length(obj.musc_obj);
            load([fileparts(pwd),'\SynergyControl\Data\Johnson2011MuscleParameters.mat'],'params');
            %We don't want to waste time going through the entire list if we've already found the muscle we're looking for.
            %Make a list of the muscle and switch off the ones we've already found
            muscle_list = ones(num_musc,1);
            for i=2:size(params,1)
                for j=1:num_musc
                    %Make sure that the Excel muscle names are in the same format as the Animatlab ones
                    name_no_space = params{i,1}{1}(~isspace(params{i,1}{1}));
                    nameAnimatlab = obj.musc_obj{j}.muscle_name;
                    if muscle_list(j,1) == 1
                        if contains(obj.musc_obj{j}.muscle_name,name_no_space)
                            %Switch this muscle off so we skip it next time
                            muscle_list(j,1) = 0;
                            %Muscle mass in [mg]
                            obj.musc_obj{j}.mass = params{i,2};
                            %Optimal muscle fiber length in [mm]
                            obj.musc_obj{j}.opt_fiber_length = params{i,3};
                            %Pinnation angle in [deg]
                            obj.musc_obj{j}.pennation_angle = params{i,4};
                            %Maximum isometric force in [g]
                            obj.musc_obj{j}.Po = params{i,5};
                            %Maximum fiber shortening velocity in [mm/s]
                            obj.musc_obj{j}.vmax_fiber = params{i,6};
                            %Tendon slack length in [mm]
                            obj.musc_obj{j}.tendon_sl = params{i,7};
                            %Optimal muscle length from Eng 2008 (cm)
                            obj.musc_obj{j}.opt_muscle_length = params{i,8};
                            %Lf/Lm, fiber length to muscle length ratio
                            obj.musc_obj{j}.lf_lm = params{i,9};
                            break
                        end
                    end
                end
            end
        end
        %% Function: Store Muscle Profiles, Attachment Points, and Muscle Length Profiles in Muscle Objects
        function store_muscle_profiles(obj)
            %Uses the joint rotation matrices at each timestep to calculate the muscle length and velocity. Stores this as a profile in the muscle object.
            num_steps = size(obj.joint_obj{1}.joint_rotmat_profile,3);
           
            attachments_body = cell(4,1);
            attdirect = attachments_body;
            attachments_post = cell(4,1,num_steps);
            
            %Increments the pointer for each body so that new attachment entries don't overwrite previous ones
            %Also used to count how many attachments each body has
            pcount = 1;
            fcount = 1;
            tcount = 1;
            ftcount = 1;
            
            %This loop builds two cell arrays: one that holds the attachment positions for each body and one that functions as a directory for storing the
            %attachment positions in the correct muscles. The directoy holds (for each body) the muscle number in the top row and the attachment number in the
            %bottom row.
            for i = 1:size(obj.musc_obj,1)
                for j = 1:size(obj.musc_obj{i}.pos_attachments,1)
                    switch obj.musc_obj{i}.pos_attachments{j,3}
                        case 1
                            attachments_body{1}(:,pcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{1}(1,pcount) = i;
                            attdirect{1}(2,pcount) = j;
                            pcount = pcount + 1;
                        case 2
                            attachments_body{2}(:,fcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{2}(1,fcount) = i;
                            attdirect{2}(2,fcount) = j;
                            fcount = fcount + 1;
                        case 3
                            attachments_body{3}(:,tcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{3}(1,tcount) = i;
                            attdirect{3}(2,tcount) = j;
                            tcount = tcount + 1;
                        case 4
                            attachments_body{4}(:,ftcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{4}(1,ftcount) = i;
                            attdirect{4}(2,ftcount) = j;
                            ftcount = ftcount + 1;
                    end
                end
            end
            
            hipJointRotMat = obj.joint_obj{1}.joint_rotmat_profile;
            kneeJointRotMat = obj.joint_obj{2}.joint_rotmat_profile;
            ankleJointRotMat = obj.joint_obj{3}.joint_rotmat_profile;
            
            orgPos = obj.organism_position';
            
            pelvisCR = obj.body_obj{1}.CR;
            femurCR = obj.body_obj{2}.CR;
            tibiaCR = obj.body_obj{3}.CR;
            footCR = obj.body_obj{4}.CR;
            
            femurPos = obj.body_obj{2}.position';
            tibiaPos = obj.body_obj{3}.position';
            footPos = obj.body_obj{4}.position';
            
            %Now to actually move the attachments through the motion. This loop builds a post-motion cell array for each body
           for i = 1:num_steps
                a = orgPos+pelvisCR*attachments_body{1};
                b = orgPos+pelvisCR*(femurPos+hipJointRotMat(:,:,i)*femurCR*attachments_body{2});
                c = orgPos+pelvisCR*(femurPos+hipJointRotMat(:,:,i)*femurCR*...
                    (tibiaPos+kneeJointRotMat(:,:,i)*tibiaCR*attachments_body{3}));
                d = orgPos+pelvisCR*(femurPos+hipJointRotMat(:,:,i)*femurCR*...
                    (tibiaPos+kneeJointRotMat(:,:,i)*tibiaCR*...
                    (footPos+ankleJointRotMat(:,:,i)*footCR*attachments_body{4})));
                attachments_post(:,:,i) = [{a};{b};{c};{d}];
            end
            
            %This loop creates an array of the attachment's motion and then stores it in the correct muscle object from the directory
            for j = 1:size(attachments_post,1)
                num_muscles = size(attachments_post{j},2);
                bb = [attachments_post{j,1,:}]';
                r = size(bb);
                for k = 1:size(attachments_post{j},2)
                    musclenumber = attdirect{j}(1,k);
                    attnumber = attdirect{j}(2,k);
                    obj.musc_obj{musclenumber}.pos_attachments{attnumber,4} = bb(k:num_muscles:r,:);
                end
            end
            
            %This loop calculates the muscle length from the attachment points and stores the length profile (and velocity) in the muscle object
            for i = 1:length(obj.musc_obj)
                %Put this muscle's attachments into a single matrix
                attachmentmatrix = cell2mat(obj.musc_obj{i}.pos_attachments(:,4));
                attatts = zeros(num_steps,size(obj.musc_obj{i}.pos_attachments,1)-1);
                for k = 1:size(obj.musc_obj{i}.pos_attachments,1)-1
                    %This just iterates through each attachment pair, subtracting the distal position from the proximal position.
                    %Since attachment positions are in one long matrix, we have to set moving pointers for the beginning and end of each attachment profile
                    temp = attachmentmatrix((k*num_steps+1):((k+1)*num_steps),:)-attachmentmatrix(((k-1)*num_steps+1):k*num_steps,:);
                    %Once we have the difference between the two attachments, take the vector norm along the second dimension (each row)
                    attatts(:,k) = vecnorm(temp,2,2);
                end
                %Each column of attatts corresponds to the length of a different segment. Summing each segment gives us the full muscle length.
                muscle_length = sum(attatts,2);
                obj.musc_obj{i}.muscle_length_profile = muscle_length;
                obj.musc_obj{i}.l_min = min(muscle_length);
                obj.musc_obj{i}.l_max = max(muscle_length);
                vTemp = diff(smoothdata(muscle_length(1:end-9),'loess',300))/obj.dt_motion;
                vTemp = [vTemp;repmat(vTemp(end),10,1)];
                obj.musc_obj{i}.muscle_velocity_profile = vTemp;
            end   
        end
        %% Function: Compute Muscle Lengths for Input Theta Profile
        function outLens = muscle_lengths_on_demand(obj,theta)
            %Uses the joint rotation matrices at each timestep to calculate the muscle length and velocity. Stores this as a profile in the muscle object.
            num_steps = size(obj.joint_obj{1}.joint_rotmat_profile,3);
           
            attachments_body = cell(4,1);
            attdirect = attachments_body;
            attachments_post = cell(4,1,num_steps);
            
            %Increments the pointer for each body so that new attachment entries don't overwrite previous ones
            %Also used to count how many attachments each body has
            pcount = 1;
            fcount = 1;
            tcount = 1;
            ftcount = 1;
            
            %This loop builds two cell arrays: one that holds the attachment positions for each body and one that functions as a directory for storing the
            %attachment positions in the correct muscles. The directoy holds (for each body) the muscle number in the top row and the attachment number in the
            %bottom row.
            for i = 1:size(obj.musc_obj,1)
                for j = 1:size(obj.musc_obj{i}.pos_attachments,1)
                    switch obj.musc_obj{i}.pos_attachments{j,3}
                        case 1
                            attachments_body{1}(:,pcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{1}(1,pcount) = i;
                            attdirect{1}(2,pcount) = j;
                            pcount = pcount + 1;
                        case 2
                            attachments_body{2}(:,fcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{2}(1,fcount) = i;
                            attdirect{2}(2,fcount) = j;
                            fcount = fcount + 1;
                        case 3
                            attachments_body{3}(:,tcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{3}(1,tcount) = i;
                            attdirect{3}(2,tcount) = j;
                            tcount = tcount + 1;
                        case 4
                            attachments_body{4}(:,ftcount) = obj.musc_obj{i}.pos_attachments{j,1};
                            attdirect{4}(1,ftcount) = i;
                            attdirect{4}(2,ftcount) = j;
                            ftcount = ftcount + 1;
                    end
                end
            end
            
            axesMat = joint_axes_from_angles(obj,theta);
            hipJointRotMat   = axis_angle_rotation(obj,theta(1),axesMat(:,1));
            kneeJointRotMat  = axis_angle_rotation(obj,theta(2),axesMat(:,2));
            ankleJointRotMat = axis_angle_rotation(obj,theta(3),axesMat(:,3));
            
            pelvisCR = obj.body_obj{1}.CR;
            femurCR = obj.body_obj{2}.CR;
            tibiaCR = obj.body_obj{3}.CR;
            footCR = obj.body_obj{4}.CR;
            
            femurPos = obj.body_obj{2}.position';
            tibiaPos = obj.body_obj{3}.position';
            footPos = obj.body_obj{4}.position';
            
            %Now to actually move the attachments through the motion. This loop builds a post-motion cell array for each body
                a = pelvisCR*attachments_body{1};
                b = pelvisCR*(femurPos+hipJointRotMat*femurCR*attachments_body{2});
                c = pelvisCR*(femurPos+hipJointRotMat*femurCR*...
                    (tibiaPos+kneeJointRotMat*tibiaCR*attachments_body{3}));
                d = pelvisCR*(femurPos+hipJointRotMat*femurCR*...
                    (tibiaPos+kneeJointRotMat*tibiaCR*...
                    (footPos+ankleJointRotMat*footCR*attachments_body{4})));
                attachments_post = [{a};{b};{c};{d}];
            
                muscCell = cell(7,1);
            %This loop creates an array of the attachment's motion and then stores it in the correct muscle object from the directory
            for j = 1:size(attachments_post,1)
                num_muscles = size(attachments_post{j},2);
                bb = [attachments_post{j,1,:}]';
                r = size(bb);
                for k = 1:size(attachments_post{j},2)
                    xx = k:num_muscles:r;
                    temp = bb(xx,:)';
                    musclenumber = attdirect{j}(1,k);
                    attnumber = attdirect{j}(2,k);
                    muscCell{musclenumber,attnumber} = temp';
                end
            end
            outLens = zeros(size(muscCell,1),1);
            %This loop calculates the muscle length from the attachment points and stores the length profile (and velocity) in the muscle object
            for i = 1:size(muscCell,1)
%                 %Put this muscle's attachments into a single matrix
%                 attachmentmatrix = cell2mat(obj.musc_obj{i}.pos_attachments(:,4));
%                 attatts = zeros(num_steps,size(obj.musc_obj{i}.pos_attachments,1)-1);
%                 for k = 1:size(obj.musc_obj{i}.pos_attachments,1)-1
%                     %This just iterates through each attachment pair, subtracting the distal position from the proximal position.
%                     %Since attachment positions are in one long matrix, we have to set moving pointers for the beginning and end of each attachment profile
%                     temp = attachmentmatrix((k*num_steps+1):((k+1)*num_steps),:)-attachmentmatrix(((k-1)*num_steps+1):k*num_steps,:);
%                     %Once we have the difference between the two attachments, take the vector norm along the second dimension (each row)
%                     attatts(:,k) = vecnorm(temp,2,2);
%                 end
                %Each column of attatts corresponds to the length of a different segment. Summing each segment gives us the full muscle length.
               temp1 = muscCell(i,:);
               temp2 = reshape(cell2mat(temp1),3,[]);
               outLens(i) = sum(vecnorm(diff(temp2,1,2),2,1));
            end   
        end
        %% Function: Store Input Muscle Passive Tension
        function store_input_muscle_passive_tension(obj,passive_tension)
            musc_names = cell(size(obj.musc_obj,1),1);
            for ii = 1:size(obj.musc_obj,1)
                musc_names{ii,1} = obj.musc_obj{ii}.muscle_name;
            end
            for ii = 1:size(passive_tension,1)
                ptName = passive_tension{ii,1};
                muscInd = find(contains(musc_names,ptName),1,'first');
                obj.musc_obj{muscInd}.passive_tension = passive_tension{ii,2};
                obj.passive_tension(:,muscInd) = passive_tension{ii,2}';
            end
        end
        %% Function: Store the sampling vector
        function store_sampling_vector(obj)
            %% Optimization and simulations can take ages if not downsampled. This function creates a sampling vector which parses down analysis to a single step
            [beg,ennd,~] = find_step_indices(obj);
            div = 500;
            obj.sampling_vector = floor(linspace(beg,ennd,div));
            
              alt = linspace(1,length(obj.theta_motion),length(obj.theta_motion));
              obj.sampling_vector = alt;
        end
        %% Function: Find the Global Point Position Profile for a Local Point (input) in a Body Coordinate Frame (input)
        function [point_profile] = move_point_through_walking(obj,bodynum,localpoint)
            %Uses the joint rotation matrices at each timestep to calculate the muscle length and velocity. Stores this as a profile in the muscle object.
            num_steps = size(obj.joint_obj{1}.joint_rotmat_profile,3);
            if size(localpoint,2) == 3
                localpoint = localpoint';
            end
            
            if log(norm(localpoint)/norm(obj.body_obj{2}.position')) > 2
                localpoint = localpoint./1000;
            end
            
            point_profile = zeros(3,num_steps);
            
            switch bodynum
                case 1
                    for i = 1:num_steps
                        point_profile(:,i) = obj.body_obj{1}.CR*localpoint;
                    end
                case 2
                    for i = 1:num_steps
                        point_profile(:,i) = obj.body_obj{1}.CR*(obj.body_obj{2}.position'+obj.joint_obj{1}.joint_rotmat_profile(:,:,i)*obj.body_obj{2}.CR*localpoint);
                    end
                case 3
                    for i = 1:num_steps
                        point_profile(:,i) = obj.body_obj{1}.CR*(obj.body_obj{2}.position'+obj.joint_obj{1}.joint_rotmat_profile(:,:,i)*obj.body_obj{2}.CR*...
                            (obj.body_obj{3}.position'+obj.joint_obj{2}.joint_rotmat_profile(:,:,i)*obj.body_obj{3}.CR*localpoint));
                    end
                case 4
                    for i = 1:num_steps
                        point_profile(:,i) = obj.body_obj{1}.CR*(obj.body_obj{2}.position'+obj.joint_obj{1}.joint_rotmat_profile(:,:,i)*obj.body_obj{2}.CR*...
                            (obj.body_obj{3}.position'+obj.joint_obj{2}.joint_rotmat_profile(:,:,i)*obj.body_obj{3}.CR*...
                            (obj.body_obj{4}.position'+obj.joint_obj{3}.joint_rotmat_profile(:,:,i)*obj.body_obj{4}.CR*localpoint)));
                    end
            end 
        end
        %% Function: Write Muscle Parameters to Animatlab
        function [parameters] = write_parameters_to_animatlab(obj)
            num_muscles = size(obj.musc_obj,1);
            parameters = cell(num_muscles,1);
            project_file = importdata(obj.proj_file);
            muscle_addresses = contains(project_file,'<Type>LinearHillMuscle</Type>');
            muscle_indices = find(muscle_addresses)-2;
            %Some paramters terms have placeholders because we want to write that parameter to our Matlab objects but not overwrite them in the simulation file
            parameter_terms = {'NamePlaceholder';...
                               '<B Value';...
                               '<Lwidth Value';...
                               'VmaxPlaceHolder';...
                               '<Kse Value';...
                               '<Kpe Value';...
                               '<B Value';...
                               '<D Value';...
                               '<LowerLimitScale Value';...
                               '<UpperLimitScale Value';...
                               '<RestingLength'};
            for i=1:num_muscles
                parameters{1,1} = 'Muscle name';
                parameters{i+1,1} = obj.musc_obj{i}.muscle_name;
                parameters{1,2} = 'Maximum force';
                parameters{i+1,2} = obj.musc_obj{i}.Po*(9.8/1000);
                parameters{1,3} = 'L_width (cm)';
                parameters{i+1,3} = calculate_muscle_width(obj,obj.musc_obj{i})*100;
                parameters{1,4} = 'V_max Muscle';
                parameters{i+1,4} = (1/obj.musc_obj{i}.lf_lm)*obj.musc_obj{i}.vmax_fiber;
                parameters{1,5} = 'Kse';
                if obj.musc_obj{i}.tendon_sl == 0
                    %This is a trendline found from comparing tendon slack length to muscle mass. There was a R^2 = .7089 correlation between the parameters.  
                    parameters{i+1,5} = 78.278*obj.musc_obj{i}.mass+6952.7;
                else
                    parameters{i+1,5} = 2.5*parameters{i+1,2}/(0.067*0.001*obj.musc_obj{i}.tendon_sl);
                end
                parameters{1,6} = 'Kpe';
                parameters{i+1,6} = (.3*parameters{i+1,5}*parameters{i+1,2})/(parameters{i+1,5}*(obj.musc_obj{i}.l_max-obj.musc_obj{i}.l_min)-.3*parameters{i+1,2});
                if parameters{i+1,6} < 0
                    disp('Something is wrong, the calculated Kpe value is negative.')
                    keyboard
                end
                %The damping constant is the maximum muscle force in N divded by the maximum velocity of the muscle in m/s
                parameters{1,7} = 'B';
                parameters{i+1,7} = (parameters{i+1,2}/parameters{i+1,4})*1000;
                parameters{1,8} = 'Yoffset (mN)';
                parameters{i+1,8} = -2.3*parameters{i+1,2};
                parameters{1,9} = 'l_min (cm)';
                parameters{i+1,9} = min(obj.musc_obj{i}.muscle_length_profile)*100;
                parameters{1,10} = 'l_max (cm)';
                parameters{i+1,10} = max(obj.musc_obj{i}.muscle_length_profile)*100;
                parameters{1,11} = 'l_rest (cm)';
                parameters{i+1,11} = parameters{i+1,10};
                for j=1:size(parameters,2)
                    lower_limit = muscle_indices(i);
                    scale = 1;
                    if j ~= 1 && j ~= 4 
                        if j == 7
                            %For some dumb reason, the creator fo Animatlab has two parameters named 'B'. In order to put damping in the right place, we have to
                            %skip over the first B.
                            prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                            lower_limit = find(prop_addresses,1,'first')+lower_limit;
                        end
                        if j == 9 || j == 10
                            %We need to do a similar skip over w Lmin and Lmax since they're used for two different things.
                            prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                            lower_limit = find(prop_addresses,1,'first')+lower_limit;
                        end
                    prop_addresses = contains(project_file(lower_limit:end),parameter_terms{j});
                    prop_index = find(prop_addresses,1,'first')+lower_limit-1;
                    %Find the line with the parameter
                    line_of_interest = project_file{prop_index};
                    if contains(line_of_interest,'None') == 1
                    else
                        if contains(line_of_interest,'milli') == 1
                            scale = 1000;
                        elseif contains(line_of_interest,'centi') == 1
                            scale = 100;
                        end
                    end
                    quotelocs = strfind(line_of_interest,'"');
                    modified_line = strcat(line_of_interest(1:quotelocs(1)),num2str(parameters{i+1,j}),line_of_interest(quotelocs(2):quotelocs(end-1)),num2str(parameters{i+1,j}/scale),line_of_interest(quotelocs(end):end));
                    %Replace the line with the modified parameter
                    project_file{prop_index} = modified_line;
                    end
                end
            end
            carry_on = input('You are about to overwrite the Animatlab project file you''re using with new parameters.\n This could permanently ruin your project file if there are errors.\n If this is what you want to do, type YES. Otherwise, the program will not overwrite.\n','s');
            if strcmp(carry_on,'YES')
                filename = 'G:\My Drive\Rat\Optimizer\CalculatedPropertiesForAnimatlab.xlsx';
                xlswrite(filename,parameters);
                %file_path = strcat(obj.proj_file(1:end-6),'_2',obj.proj_file(end-5:end));
                file_path = strcat(obj.proj_file);
                fileID = fopen(file_path,'w');
                formatSpec = '%s\n';
                nrows = size(project_file);
                for row = 1:nrows
                    fprintf(fileID,formatSpec,project_file{row,:});
                end
                fclose(fileID);
            end
        end
        %% Function: Calculate the rotation matrix for an arbitrary set of three angles
        function C_mat = three_axis_rotation(obj,angs)
            %Take in a euler angle triplet, and return the associated
            %rotation matrix.
            c1 = cos(angs(1));
            s1 = sin(angs(1));
            c2 = cos(angs(2));
            s2 = sin(angs(2));
            c3 = cos(angs(3));
            s3 = sin(angs(3));
            %XYZ            
            C_mat = [c2*c3,-c2*s3,s2;...
                     c1*s3+c3*s1*s2,c1*c3-s1*s2*s3,-c2*s1;...
                     s1*s3-c1*c3*s2,c3*s1+c1*s2*s3,c1*c2];
        end
        %% Function: Calculate Rotation Matrix for an Input Rotation Angle
        function C = axis_angle_rotation(obj,angle,joint_axis)
            c = cos(angle);
            s = sin(angle);

            a1 = joint_axis(1);
            a2 = joint_axis(2);
            a3 = joint_axis(3);
            
            C = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
                 a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
                 a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];
        end
        %% Function: FIND find_step_indices: STEP INDICES Find the bounding indices for a single step
        function [beg,ennd,mid] = find_step_indices(obj)
            %Finds the bounding indices for the second step in walking. This will be used to find muscle moment arms and passive tension
            %Over a step since individual joint motion can't be isolated. Using the second step prevents issues with initialization of the first step (t==0)

            jointprofile = obj.theta_motion;
            trimmedjp = jointprofile(floor(length(jointprofile)*.1):floor(length(jointprofile)*.9),:);
            beg = 1;
            ennd = size(jointprofile,1);
            mid = floor((beg+ennd)/2);
            rangeJP = max(trimmedjp,[],'all')-min(trimmedjp,[],'all');
            if mean(rangeJP) < 1e-3 || mean(std(trimmedjp)) < .01
                % If jointprofile is constant
                beg = floor(length(jointprofile)*(1/3));
                ennd = floor(length(jointprofile)*(2/3));
                mid = floor((beg+ennd)/2);
                return
            end
            
            if size(unique(jointprofile(:,1)),1) == 1
                if size(unique(jointprofile(:,2)),1) == 1
                    maxvals = findpeaks(jointprofile(:,3));
                    [minvals,minlocs]= findpeaks(-jointprofile(:,3));
                else
                    maxvals = findpeaks(jointprofile(:,2));
                    [minvals,minlocs]= findpeaks(-jointprofile(:,2));
                end
            else
                maxvals = findpeaks(trimmedjp(:,1));
                [minvals,minlocs]= findpeaks(-trimmedjp(:,1));
            end
            if size(maxvals,1) > 2
%                 [minvals,minlocs]= findpeaks(-jointprofile(:,1));
                minvals = -minvals;
                if size(minvals,1) == size(maxvals,1)
                    beg = minlocs(1);
                    ennd = minlocs(2);
                else
                    beg = minlocs(end-1);
                    ennd = minlocs(end);
                end
                mid = floor((beg+ennd)/2);
            elseif size(maxvals,1) == 1 && (size(minvals,1) == 1 || isempty(minvals))
                beg = 1;
                ennd = size(jointprofile,1);
                mid = floor((beg+ennd)/2);
            else
                %plot(jointprofile)
                %disp('Are you using a simulation that has at least three steps? Check that the boundaries are going to be defined correctly for a single step')
                return
            end
        end
        %% Function: MUSCLE MOMENT ARM Compute a single muscle's moment arm over an entire walking cycle
        function [moment_arm_profile] = compute_muscle_moment_arm(obj,muscle,axis_used,joint,plotval)
            if isnumeric(muscle)
                muscle = obj.musc_obj{muscle};
            end
            
            if muscle.pos_attachments{1,3} > joint || muscle.pos_attachments{end,3} < joint
                disp('This muscle may not bridge the joint specified')
                moment_arm_profile = -1;
                return
            end
            
            if nargin < 4 || ~(plotval == 1 || plotval == 0)
                plotval = 0;
            end

            scale = 1000;
            
            if joint >= 1 && joint <= 3
                jointObj = obj.joint_obj{joint};
            else
                joint = input('Joint specified incorrectly. hip = 1, knee = 2, ankle = 3');
            end

%             [beg,ennd,~] = find_step_indices(obj);
%             %Higher div, the longer the program takes.
%             div = 500;
%             plotX = floor(linspace(beg,ennd,div));
%             xx = linspace(1,length(obj.theta_motion),length(obj.theta_motion));
            xx = obj.sampling_vector;
            %%% REVERT THIS OR THE PROGRAM WILL TAKE AGES
            %xx = 1:length(obj.theta_motion);

            axis = ['Ext/Flx';'Abd/Add';'ExR/InR'];
            %Defining values for plot axis limits
            if plotval == 1 
                jointmat = [obj.joint_obj{1}.sim_position_profile(xx,:);...
                                obj.joint_obj{2}.sim_position_profile(xx,:);...
                                obj.joint_obj{3}.sim_position_profile(xx,:)];
                if size(obj.musc_obj,1) ~= 38
                    toeMusc = obj.musc_obj{6,1};
                else
                    toeMusc = obj.musc_obj{20,1};
                end
                toemat = toeMusc.pos_attachments{5,4};
                jointaxmat = (1/100).*[jointObj.uuw_joint(:,xx)';...
                    jointObj.uuw_joint2(:,xx)';...
                    jointObj.uuw_joint3(:,xx)'];
                
                muscmax = max(cell2mat(muscle.pos_attachments(:,4)));
                muscmin = min(cell2mat(muscle.pos_attachments(:,4)));
                jointmax = max(jointmat);
                jointmin = min(jointmat);
                toemax = max(toemat);
                toemin = min(toemat);
                jointaxmax = max(jointaxmat);
                jointaxmin = min(jointaxmat);

                posmax = scale*max([jointmax;...
                                     muscmax;...
                                     toemax;...
                                     jointaxmax]);
                posmin = scale*min([jointmin;...
                                    muscmin;...
                                    toemin;...
                                    jointaxmin]);
                limits =  [posmin' posmax'];
                xlims = limits(1,:);
                ylims = limits(2,:);
                zlims = limits(3,:);
                %Can change limscale to give more or less viewing space
                limscale = 1;
            end 
                moment_arm_profile = zeros(size(xx,2),5);
                moment_arm_profile(:,1) = jointObj.rec_angle_time(xx);
                moment_arm_profile(:,2) = jointObj.rec_angle_profile(xx)*(180/pi);

            for j = axis_used
                %f1 = figure(1);
                %for i=1:length(jointprofile)
                count2 = 0;
                whichsegmentisfree = diff(cell2mat(muscle.pos_attachments(:,3)));
                if size(find(whichsegmentisfree>0),1) > 1
                    atts = cell2mat(muscle.pos_attachments(:,3));
                    holder  = atts-joint;
                    seg2use = find(holder==1,1,'first')-1;
                    if whichsegmentisfree(seg2use,1) == 0
                        keyboard
                    end
                else
                    [~,seg2use] = max(whichsegmentisfree);
                end
                jointvBook = (scale/100)*[reshape(jointObj.uuw_joint,3,1,size(jointObj.uuw_joint,2)),...
                                          reshape(jointObj.uuw_joint2,3,1,size(jointObj.uuw_joint2,2)),...
                                          reshape(jointObj.uuw_joint3,3,1,size(jointObj.uuw_joint3,2))];
                a = muscle.pos_attachments{seg2use+1-1,4} - jointObj.sim_position_profile;
                b = muscle.pos_attachments{seg2use+2-1,4} - jointObj.sim_position_profile;
                %pointrelBook = reshape([reshape(a,1,length(a)*3);reshape(b,1,length(b)*3)],length(a)*2,3);
                switch axis_used
                    case 1
                        ppjoint = (scale/100)*jointObj.uuw_joint;
                    case 2
                        ppjoint = (scale/100)*jointObj.uuw_joint2;
                    case 3
                        ppjoint = (scale/100)*jointObj.uuw_joint3;
                end
                c = scale.*(a-(dot(a',ppjoint)'./vecnorm(ppjoint).^2').*ppjoint')+scale*jointObj.sim_position_profile;
                d = scale.*(b-(dot(b',ppjoint)'./vecnorm(ppjoint).^2').*ppjoint')+scale*jointObj.sim_position_profile;
                pointprojectionBook = reshape([reshape(c,1,length(c)*3);reshape(d,1,length(d)*3)],length(c)*2,3);
                %Matrix of muscle projections perpendicular to the joint axis
                fflatmatBook = d-c;
                momentarmlongBook = cross(fflatmatBook,ppjoint');
                scaledJointBook = scale*jointObj.sim_position_profile;
                PA2Book = reshape([reshape(c,1,length(c)*3);reshape(scaledJointBook,1,length(scaledJointBook)*3)],length(c)*2,3);
                PB2Book = reshape([reshape(d,1,length(d)*3);reshape(scaledJointBook+momentarmlongBook,1,length(scaledJointBook+momentarmlongBook)*3)],length(d)*2,3);
                moment_armBook = lineIntersect3D(obj,PA2Book,PB2Book);
                sig_momentarm = moment_armBook-scaledJointBook;
                sig_muscleprojection = c-d;
                signal2 = sign(dot(cross(sig_momentarm,sig_muscleprojection),squeeze(jointvBook(:,j,:))',2));
                moment_arm_profile(:,j+2) = signal2.*vecnorm(sig_momentarm,2,2);
                if plotval == 1
                    %% For plotting the muscle and joint in 3D as the muscle moves
                    for i = xx 
                        count2 = count2 + 1;
    %                         sig_momentarmb = moment_armBook(i,:)-scaledJointBook(i,:);
    %                         sig_muscleprojectionb = pointprojectionBook(2*i-1,:)-pointprojectionBook(2*i,:);
    %                         signal2b = sign(dot(cross(sig_momentarmb,sig_muscleprojectionb),jointvBook(:,j,i)'));
    %                         %moment_arm_length2(i-beg+1,j+1) = momentlen2;
    %                         moment_arm_profileb(count2,j+2) = signal2b*norm(sig_momentarmb);
                        %Unnecessary to plot every time step. This sets it to plot only after a certain number of steps
                        fps = 32;
                        if plotval == 1 && any(count2==1:fps:length(xx))
                            if i == 1
                                figure('name','legplot','Position',[2,2,958,994]);
                            end
                            %The muscle vector
                            for k = 1:size(muscle.pos_attachments,1)
                                musclevecco(k,:) = scale*muscle.pos_attachments{k,4}(i,:);
                            end
                            jointvecco1 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,1,i)'];
                            jointvecco2 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,2,i)'];
                            jointvecco3 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,3,i)'];
                            femur = scale*[obj.joint_obj{1}.sim_position_profile(i,1),obj.joint_obj{1}.sim_position_profile(i,2),obj.joint_obj{1}.sim_position_profile(i,3);...
                                obj.joint_obj{2}.sim_position_profile(i,1),obj.joint_obj{2}.sim_position_profile(i,2),obj.joint_obj{2}.sim_position_profile(i,3)];
                            tibia = scale*[obj.joint_obj{2}.sim_position_profile(i,1),obj.joint_obj{2}.sim_position_profile(i,2),obj.joint_obj{2}.sim_position_profile(i,3);...
                                obj.joint_obj{3}.sim_position_profile(i,1),obj.joint_obj{3}.sim_position_profile(i,2),obj.joint_obj{3}.sim_position_profile(i,3)];
                            foot = scale*[obj.joint_obj{3}.sim_position_profile(i,1),obj.joint_obj{3}.sim_position_profile(i,2),obj.joint_obj{3}.sim_position_profile(i,3);...
                                toemat(i,:)];
                            %Plot the muscle
                            gcf;
                            plot3(musclevecco(:,1),musclevecco(:,2),musclevecco(:,3),'r', 'LineWidth', 2)
                            %Thicken the muscle
    %                         set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', 2);
                            hold on
                            %Plot the moment arm average
                            singlesegmentmomentarm = [scaledJointBook(i,:);moment_armBook(i,:)];
                            plot3(singlesegmentmomentarm(:,1),singlesegmentmomentarm(:,2),singlesegmentmomentarm(:,3),'m','LineWidth',4.5) 

                            % Plot the long moment arm
    %                         longmomentarm = [scaledJointBook(i,:);momentarmlongBook(i,:)];
    %                         plot3(longmomentarm(:,1),longmomentarm(:,2),longmomentarm(:,3),'m','LineWidth',4.5)

                            %Plot the muscle origin
                            plot3(scale*muscle.pos_attachments{1,4}(i,1),scale*muscle.pos_attachments{1,4}(i,2),scale*muscle.pos_attachments{1,4}(i,3),'ro')
                            %Plot the muscle insertion
                            plot3(scale*muscle.pos_attachments{end,4}(i,1),scale*muscle.pos_attachments{end,4}(i,2),scale*muscle.pos_attachments{end,4}(i,3),'r+')
                            %Plot the joint of interest
                            plot3(scale*obj.joint_obj{joint}.sim_position_profile(i,1),scale*obj.joint_obj{joint}.sim_position_profile(i,2),scale*obj.joint_obj{joint}.sim_position_profile(i,3),'b.')
                            %Plot the Ext/Flx axis
                            plot3(jointvecco1(:,1),jointvecco1(:,2),jointvecco1(:,3),'r')
                            %Plot the Add/Abd axis
                            plot3(jointvecco2(:,1),jointvecco2(:,2),jointvecco2(:,3),'g')
                            %Plot the ExR/InR axis
                            plot3(jointvecco3(:,1),jointvecco3(:,2),jointvecco3(:,3),'b')
                            %Plot the muscle projection on the plane of the axis of interest
                            plot3(pointprojectionBook(2*i-1:2*i,1),pointprojectionBook(2*i-1:2*i,2),pointprojectionBook(2*i-1:2*i,3),'m--')
                            %Plot the hip
                            jointSize = 12;
                            plot3(scale*obj.joint_obj{1}.sim_position_profile(i,1),scale*obj.joint_obj{1}.sim_position_profile(i,2),scale*obj.joint_obj{1}.sim_position_profile(i,3),'kp','MarkerSize',jointSize)
                            %Plot the femur
                            boneWidth = 3;
                            plot3(femur(:,1),femur(:,2),femur(:,3),'k--','LineWidth',boneWidth)
                            %Plot the knee
                            plot3(scale*obj.joint_obj{2}.sim_position_profile(i,1),scale*obj.joint_obj{2}.sim_position_profile(i,2),scale*obj.joint_obj{2}.sim_position_profile(i,3),'ks','MarkerSize',jointSize)
                            %Plot the tibia
                            plot3(tibia(:,1),tibia(:,2),tibia(:,3),'k-.','LineWidth',boneWidth)
                            %Plot the ankle
                            plot3(scale*obj.joint_obj{3}.sim_position_profile(i,1),scale*obj.joint_obj{3}.sim_position_profile(i,2),scale*obj.joint_obj{3}.sim_position_profile(i,3),'kd','MarkerSize',jointSize)
                            %Plot the foot
                            plot3(foot(:,1),foot(:,2),foot(:,3),'k:','LineWidth',boneWidth)
                            hold on
                            grid on
                            %Set figure properties
                            title([{[muscle.muscle_name(4:end),' ',axis(j,:),' about the ',obj.joint_obj{joint}.name(4:end),' joint']},...
                                   {['Moment arm length: ',num2str(round(moment_arm_profile(count2,j+2),2)),' mm']},...
                                   {[obj.joint_obj{joint}.name(4:end),' joint angle: ',num2str(obj.joint_obj{joint}.rec_angle_profile(i)*(180/pi))]}])
                            legend('Muscle','Moment arm','Muscle Origin','Muscle Insertion','Joint of Interest','Joint vector Ext/Flx','Joint vector Add/Abd','Joint vector ExR/InR',...
                                'Muscle Projection onto Joint Plane','Hip','Femur','Knee','Tibia','Ankle','Foot')
                            view([0 90])
    %                         xlim([-20 35])
    %                         ylim([-60 10])
                            xlim(limscale.*xlims)
                            ylim(limscale.*ylims)
                            zlim([zlims(1)-sign(zlims(1))*limscale*zlims(1) zlims(2)+sign(zlims(2))*limscale*zlims(2)])
                            xlabel('X')
                            ylabel('Y')
                            zlabel('Z')
                            %To avoid distorted 3D plots, set the aspect ratio of the plot relative to the lengths of the axes
                            axeslengths = [range(xlim);range(ylim);range(zlim)];
                            normedaxes = axeslengths/norm(axeslengths);
                            pbaspect(normedaxes)
                            hold off
                            pause(.05)
                            if i == fps+1
                                pause
                            end
                        end
                    end 
                end
            end
        end
        %% Function: GRAVITY MOMENT ARM: Compute the gravity moment arm for a given body or joint
        function [moment_arm_profile] = compute_gravity_moment_arm(obj,inObj,joint,plotval)
            
            objName = inObj.name;
            body_names = cell(1,length(obj.body_obj));
            joint_names = cell(1,length(obj.joint_obj));
            for ii = 1:length(obj.body_obj)
                body_names{ii} = obj.body_obj{ii}.name;
                if ii < length(obj.body_obj)
                    joint_names{ii} = obj.joint_obj{ii}.name;
                end
            end
            if any(ismember(joint_names,objName))
                jointInd = find(contains(joint_names,objName),1,'first');
                if isempty(jointInd)
                    error('compute_gravity_moment_arm: Joint object not found in FullLeg.')
                elseif joint == jointInd
                    error('compute_gravity_moment_arm: Provided joint will not generate a moment arm about itself.')
                else
                    point_profile = obj.joint_obj{jointInd}.sim_position_profile';
                end
            elseif any(ismember(body_names,objName))
                bodyInd = find(contains(body_names,objName),1,'first');
                if isempty(bodyInd)
                    error('compute_gravity_moment_arm: Body object not found in FullLeg.')
                else
                    point_profile = move_point_through_walking(obj,bodyInd,obj.body_obj{bodyInd}.com);
                end
            else 
                error('compute_gravity_moment_arm: inObj must be a body or a joint.')
            end
           
            antiGrav = point_profile-[0 -.005 0]';
            
            if nargin < 4 || ~ismember(plotval,[0 1])
                plotval = 0;
            end

            scale = 1000;
            
            if joint >= 1 && joint <= 3
                jointObj = obj.joint_obj{joint};
            else
                joint = input('Joint specified incorrectly. hip = 1, knee = 2, ankle = 3');
            end

            [beg,ennd,~] = find_step_indices(obj);
            xx = obj.sampling_vector;

            axis = ['Ext/Flx';'Abd/Add';'ExR/InR'];
            %Defining values for plot axis limits[moment_arm_profile] = compute_body_moment_arm(obj,obj.body_obj{2},1,1,0);
            if plotval == 1 
                jointmat = [obj.joint_obj{1}.sim_position_profile(beg:ennd,:);...
                                obj.joint_obj{2}.sim_position_profile(beg:ennd,:);...
                                obj.joint_obj{3}.sim_position_profile(beg:ennd,:)];
                toemat = obj.musc_obj{20,1}.pos_attachments{5,4};
                jointaxmat = (1/100).*[obj.joint_obj{joint}.uuw_joint(:,xx)';...
                    obj.joint_obj{joint}.uuw_joint2(:,xx)';...
                    obj.joint_obj{joint}.uuw_joint3(:,xx)'];
                
                pointsmax = max(point_profile,[],2)';
                pointsmin = min(point_profile,[],2)';
                jointmax = max(jointmat);
                jointmin = min(jointmat);
                toemax = max(toemat);
                toemin = min(toemat);
                jointaxmax = max(jointaxmat);
                jointaxmin = min(jointaxmat);

                posmax = scale*max([jointmax;...
                                     pointsmax;...
                                     toemax;...
                                     jointaxmax]);
                posmin = scale*min([jointmin;...
                                    pointsmin;...
                                    toemin;...
                                    jointaxmin]);
                limits =  [posmin' posmax'];
                xlims = limits(1,:);
                ylims = limits(2,:);
                zlims = limits(3,:);
                %Can change limscale to give more or less viewing space
                limscale = 1;
            end 
                moment_arm_profile = zeros(size(xx,2),5);
                moment_arm_profile(:,1) = obj.joint_obj{joint}.rec_angle_time(xx);
                moment_arm_profile(:,2) = obj.joint_obj{joint}.rec_angle_profile(xx)*(180/pi);
                jointShifts = [98.4373 102.226 116.2473];
                jointAngleProfile = (180/pi).*obj.joint_obj{joint}.rec_angle_profile+jointShifts(joint);

            axis_used = 1;
            for j = axis_used
                count2 = 0;
                jointvBook = (scale/100)*[reshape(jointObj.uuw_joint,3,1,size(jointObj.uuw_joint,2)),...
                                          reshape(jointObj.uuw_joint2,3,1,size(jointObj.uuw_joint2,2)),...
                                          reshape(jointObj.uuw_joint3,3,1,size(jointObj.uuw_joint3,2))];
                a = antiGrav' - jointObj.sim_position_profile;
                b = point_profile' - jointObj.sim_position_profile;
                pointrelBook = reshape([reshape(a,1,length(a)*3);reshape(b,1,length(b)*3)],length(a)*2,3);
                switch axis_used
                    case 1
                        ppjoint = (scale/100)*jointObj.uuw_joint;
                    case 2
                        ppjoint = (scale/100)*jointObj.uuw_joint2;
                    case 3
                        ppjoint = (scale/100)*jointObj.uuw_joint3;
                end
                c = scale.*(a-(dot(a',ppjoint)'./vecnorm(ppjoint).^2').*ppjoint')+scale*jointObj.sim_position_profile;
                d = scale.*(b-(dot(b',ppjoint)'./vecnorm(ppjoint).^2').*ppjoint')+scale*jointObj.sim_position_profile;
                pointprojectionBook = reshape([reshape(c,1,length(c)*3);reshape(d,1,length(d)*3)],length(c)*2,3);
                %Matrix of muscle projections perpendicular to the joint axis
                fflatmatBook = d-c;
                momentarmlongBook = cross(fflatmatBook,ppjoint');
                scaledJointBook = scale*jointObj.sim_position_profile;
                PA2Book = reshape([reshape(c,1,length(c)*3);reshape(scaledJointBook,1,length(scaledJointBook)*3)],length(c)*2,3);
                PB2Book = reshape([reshape(d,1,length(d)*3);reshape(scaledJointBook+momentarmlongBook,1,length(scaledJointBook+momentarmlongBook)*3)],length(d)*2,3);
                moment_armBook = lineIntersect3D(obj,PA2Book,PB2Book);
                for i=xx
                    count2 = count2 + 1;
                        %Storing the moment arms in a matrix
                        %A modifier to determine whether the moment arm is positive or negative
                        sig_momentarm = moment_armBook(i,:)-scaledJointBook(i,:);
                        sig_muscleprojection = pointprojectionBook(2*i-1,:)-pointprojectionBook(2*i,:);
                        signal2 = sign(dot(cross(sig_momentarm,sig_muscleprojection),jointvBook(:,j,i)'));
                        %moment_arm_length2(i-beg+1,j+1) = momentlen2;
                        moment_arm_profile(count2,j+2) = signal2*norm(sig_momentarm);
                    %% For plotting the muscle and joint in 3D as the muscle moves
                    %Unnecessary to plot every time step. This sets it to plot only after a certain number of steps
                    % Higher is faster
                    fps = 200;
                    if plotval == 1 && any(count2==1:fps:length(xx))
                        if i == 1
                            figure('name','legplot','Position',[2,2,958,994]);
                        end
                        %The muscle vector
                        for k = 1:2
                            pointvecco(k,:) = scale*point_profile(:,i);
                        end
                        jointvecco1 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,1,i)'];
                        jointvecco2 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,2,i)'];
                        jointvecco3 = [scaledJointBook(i,:);scaledJointBook(i,:)+jointvBook(:,3,i)'];
                        femur = scale*[obj.joint_obj{1}.sim_position_profile(i,1),obj.joint_obj{1}.sim_position_profile(i,2),obj.joint_obj{1}.sim_position_profile(i,3);...
                            obj.joint_obj{2}.sim_position_profile(i,1),obj.joint_obj{2}.sim_position_profile(i,2),obj.joint_obj{2}.sim_position_profile(i,3)];
                        tibia = scale*[obj.joint_obj{2}.sim_position_profile(i,1),obj.joint_obj{2}.sim_position_profile(i,2),obj.joint_obj{2}.sim_position_profile(i,3);...
                            obj.joint_obj{3}.sim_position_profile(i,1),obj.joint_obj{3}.sim_position_profile(i,2),obj.joint_obj{3}.sim_position_profile(i,3)];
                        foot = scale*[obj.joint_obj{3}.sim_position_profile(i,1),obj.joint_obj{3}.sim_position_profile(i,2),obj.joint_obj{3}.sim_position_profile(i,3);...
                            obj.musc_obj{20, 1}.pos_attachments{5,4}(i,:)];
                                %flongco = [fflat(:,i-beg+1)'+scale*muscle.pos_attachments{1,4}(i,:);flong(:,i-beg+1)'+fflat(:,i-beg+1)'+scale*muscle.pos_attachments{1,4}(i,:)];
                                %fflatco = [scale*muscle.pos_attachments{1,4}(i,:);fflat(:,i-beg+1)'+scale*muscle.pos_attachments{1,4}(i,:)];
                                %mprojection = [originprojection';insprojection'];
                        %Plot the muscle
                        plot3(pointvecco(:,1),pointvecco(:,2),pointvecco(:,3),'r', 'LineWidth', 2)
                        hold on
                        %Plot the moment arm average
                        singlesegmentmomentarm = [scaledJointBook(i,:);moment_armBook(i,:)];
                        plot3(singlesegmentmomentarm(:,1),singlesegmentmomentarm(:,2),singlesegmentmomentarm(:,3),'m','LineWidth',2)                  
                        
                        %Plot the muscle origin
                        plot3(scale*antiGrav(1,i),scale*antiGrav(2,i),scale*antiGrav(3,i),'ro')
                        %Plot the muscle insertion
                        plot3(scale*point_profile(1,i),scale*point_profile(2,i),scale*point_profile(3,i),'r+')
                        %Plot the joint of interest
                        plot3(scale*obj.joint_obj{joint}.sim_position_profile(i,1),scale*obj.joint_obj{joint}.sim_position_profile(i,2),scale*obj.joint_obj{joint}.sim_position_profile(i,3),'b.')
                        %Plot the Ext/Flx axis
                        plot3(jointvecco1(:,1),jointvecco1(:,2),jointvecco1(:,3),'r')
                        %Plot the Add/Abd axis
                        plot3(jointvecco2(:,1),jointvecco2(:,2),jointvecco2(:,3),'g')
                        %Plot the ExR/InR axis
                        plot3(jointvecco3(:,1),jointvecco3(:,2),jointvecco3(:,3),'b')
                        %Plot the muscle projection on the plane of the axis of interest
                        plot3(pointprojectionBook(2*i-1:2*i,1),pointprojectionBook(2*i-1:2*i,2),pointprojectionBook(2*i-1:2*i,3),'m--')
                        %Plot the hip
                        plot3(scale*obj.joint_obj{1}.sim_position_profile(i,1),scale*obj.joint_obj{1}.sim_position_profile(i,2),scale*obj.joint_obj{1}.sim_position_profile(i,3),'kp')
                        %Plot the femur
                        plot3(femur(:,1),femur(:,2),femur(:,3),'k--')
                        %Plot the knee
                        plot3(scale*obj.joint_obj{2}.sim_position_profile(i,1),scale*obj.joint_obj{2}.sim_position_profile(i,2),scale*obj.joint_obj{2}.sim_position_profile(i,3),'ks')
                        %Plot the tibia
                        plot3(tibia(:,1),tibia(:,2),tibia(:,3),'k-.')
                        %Plot the ankle
                        plot3(scale*obj.joint_obj{3}.sim_position_profile(i,1),scale*obj.joint_obj{3}.sim_position_profile(i,2),scale*obj.joint_obj{3}.sim_position_profile(i,3),'kd')
                        %Plot the foot
                        plot3(foot(:,1),foot(:,2),foot(:,3),'k:')
                        hold on
                        grid on
                        %Set figure properties
                        title([{[objName,' ',axis(j,:),' about the ',obj.joint_obj{joint}.name(4:end),' joint']},...
                               {['Moment arm length: ',num2str(round(moment_arm_profile(count2,j+2),2)),' mm']},...
                               {[obj.joint_obj{joint}.name(4:end),' joint angle: ',num2str(jointAngleProfile(i))]}],'Interpreter','None')
                        legend('Muscle','Moment arm','Muscle Origin','Muscle Insertion','Joint of Interest','Joint vector Ext/Flx','Joint vector Add/Abd','Joint vector ExR/InR',...
                            'Muscle Projection onto Joint Plane','Hip','Femur','Knee','Tibia','Ankle','Foot')
                        view([0 90])
%                         xlim([-20 35])
%                         ylim([-60 10])
                        xlim(limscale.*xlims)
                        ylim(limscale.*ylims)
                        zlim([zlims(1)-sign(zlims(1))*limscale*zlims(1) zlims(2)+sign(zlims(2))*limscale*zlims(2)])
                        xlabel('X')
                        ylabel('Y')
                        zlabel('Z')
                        %To avoid distorted 3D plots, set the aspect ratio of the plot relative to the lengths of the axes
                        axeslengths = [range(xlim);range(ylim);range(zlim)];
                        normedaxes = axeslengths/norm(axeslengths);
                        pbaspect(normedaxes)
                        hold off
                        pause(.05)
                    end
                end 
            end
        end
        %% Function: JOINT MOMENT ARM: Compute the relevant moment arms of all muscles articulating a designated joint for a specified configuration
        function [moment_output] = compute_joint_moment_arms(obj,joint,axis)
            
            if sum([joint == 1,joint == 2,joint == 3]) ~= 1
                error('Joint must be 1, 2, or 3')
            end
            relevant_muscles = [];
            num_muscles = length(obj.musc_obj);
            
            for i = 1:num_muscles
                attachment_bodies = cell2mat(obj.musc_obj{i}.pos_attachments(:,3));
                switch joint
                    case 1
                        if attachment_bodies(1) == 1
                            relevant_muscles = [relevant_muscles;i];
                        end
                    case 2
                        if attachment_bodies(end) == 3 || (attachment_bodies(end)==4 && attachment_bodies(1)==2)
                            relevant_muscles = [relevant_muscles;i];
                        end
                    case 3
                        if attachment_bodies(end) == 4
                            relevant_muscles = [relevant_muscles;i];
                        end
                end
            end
            
            for i=1:length(relevant_muscles)
                muscle_num = relevant_muscles(i);
                moment_arm_profile = compute_muscle_moment_arm(obj,obj.musc_obj{muscle_num},axis,joint,0);
                moment_output(muscle_num,:) = moment_arm_profile(:,axis+2)';
            end
            
            moment_output(num_muscles+1,:) = moment_arm_profile(:,1)';
            moment_output(num_muscles+2,:) = moment_arm_profile(:,2)';
            
            plotter = 0;
            if plotter
                figure
                count = 1;
                CM = hsv(length(relevant_muscles));
                legendcell = {};
                for i=1:num_muscles
                    if moment_output(i,1) ~=0
                        plot(moment_output(i,:),'color',CM(count,:),'LineWidth',1.5)
                        legendcell{end+1} = obj.musc_obj{i}.muscle_name(4:end);
                        count = count +1;
                        hold on
                    end
                end
                    title(['Muscle moment arms about the ',obj.joint_obj{joint}.name(4:end)])
                    xlabel('Percent Stride')
                    ylabel('Moment arm length (mm)')
                    legend(legendcell,'Interpreter','none','Location','eastoutside')
                    set(gcf,'Position',[500 500 900 500])
                    saveas(gcf,['G:\My Drive\Rat\Optimizer\OutputFigures\moment_arms','\muscle_moment_arms','_',obj.joint_obj{joint}.name(4:end),'_',datestr(datetime('now'),'yyyymmdd'),'.png'])
            end
        end
        %% Function: On Demand: Calculate Leg Moment Arms for Given Joint Angle Vector
        function momentMat = leg_moment_arms(obj,theta)

            scale = 1000;

            numMuscles = length(obj.musc_obj);
            relevantMuscles = zeros(numMuscles,3);
            momentMat = relevantMuscles;
            
            % Create a boolean array that determines whether a muscle spans a specific joint
            % Row is muscle number, column is joint (HKA)
            for i = 1:numMuscles
                attachmentBodies = cell2mat(obj.musc_obj{i}.pos_attachments(:,3));
                relevantMuscles(i,:) = [attachmentBodies(1) == 1,attachmentBodies(end) == 3 || (attachmentBodies(end)==4 && attachmentBodies(1)==2),attachmentBodies(end) == 4];
            end

            jointv = joint_axes_from_angles(obj,theta)*10;
            jointMat = joint_pos_on_demand(obj,theta);
            
            for k = 1:numMuscles
                muscle = obj.musc_obj{k};
                for jointNum = 1:3
                    if relevantMuscles(k,jointNum) == 1
                        whichsegmentisfree = diff(cell2mat(muscle.pos_attachments(:,3)));
                        if size(find(whichsegmentisfree>0),1) > 1
                            atts = cell2mat(muscle.pos_attachments(:,3));
                            holder  = atts-jointNum;
                            seg2use = find(holder==1,1,'first')-1;
                            if whichsegmentisfree(seg2use,1) == 0
                                keyboard
                            end
                        else
                            [~,seg2use] = max(whichsegmentisfree);
                        end
                        pointRel = zeros(2,3);
                        pointProjection = pointRel;
                        pointRel = [att_pos_on_demand(obj,theta,muscle.pos_attachments(seg2use,:)) att_pos_on_demand(obj,theta,muscle.pos_attachments(seg2use+1,:))]'-jointMat(:,jointNum)';
                        pointProjection(1,:) = (pointRel(1,:) - (dot(pointRel(1,:),jointv(:,jointNum)')/norm(jointv(:,jointNum)')^2)*jointv(:,jointNum)')'+jointMat(:,jointNum);
                        pointProjection(2,:) = (pointRel(2,:) - (dot(pointRel(2,:),jointv(:,jointNum)')/norm(jointv(:,jointNum)')^2)*jointv(:,jointNum)')'+jointMat(:,jointNum);
                        pointProjection = scale*pointProjection;
                        %Matrix of muscle projections perpendicular to the joint axis
                            fflatmat = pointProjection(2,:)-pointProjection(1,:);
                        %Storing the moment arms in a matrix
                            momentArmLong = cross(fflatmat,jointv(:,jointNum));
                            scaledjoint = scale*jointMat(:,jointNum)';
                            PA2 = [pointProjection(1,:);scaledjoint];
                            PB2 = [pointProjection(2,:);scaledjoint+momentArmLong];
                        %lineINtersect3D gives the scaled vector of the moment arms for each segment
                            momentArm = lineIntersect3D(obj,PA2,PB2);
                        %A modifier to determine whether the moment arm is positive or negative
                            sig_momentarm = momentArm-scaledjoint;
                            sig_muscleprojection = (pointProjection(1,:)-pointProjection(2,:));
                            signal2 = sign(dot(cross(sig_momentarm,sig_muscleprojection),jointv(:,jointNum)'));
                            momentMat(k,jointNum) = signal2*norm(sig_momentarm);
                    end
                end
            end
        end
        %% Function: MUSCLE PASSIVE TENSION Compute a single muscle's passive tension over a joint's motion
        function [T_find,T_out,T] = compute_muscle_passive_tension(obj,muscle,length_find)
            
            %Higher div, the longer the program takes.
            div = 100;
            count = 0;
            
            %Find the maximum change in muscle length (from lowest trough to highest peak)
             muscleprofile = muscle.muscle_length_profile;
%             [~,maxlocs] = findpeaks(muscleprofile);
%             [~,minlocs] = findpeaks(-muscleprofile);
%             %minvals = minvals*-1;
%             a3 = unique([maxlocs;minlocs]);
%             a4 = muscleprofile(a3);
%             for i = 1:length(a4)-1
%                 %signer(i,1) = sign(a4(i+1,1)-a4(i,1));
%                 range(i,1) = abs(a4(i+1,1)-a4(i,1));
%             end
%             [~,dd] = max(range);
%             dd = dd(1);
%             beg = a3(dd);
%             ennd = a3(dd+1);
            [beg,ennd,~] = obj.find_step_indices;
            xx = floor(linspace(beg,ennd,div));
            xx = linspace(1,length(obj.theta_motion),length(obj.theta_motion));
            xx = obj.sampling_vector;
            
            T = zeros(length(xx),1);
            delx = T;
            Tdot = T;
            mV = T;
            
            ks = muscle.Kse;
            kp = muscle.Kpe;
            dt = obj.dt_motion;
            b = muscle.damping;
            Lr = muscle.RestingLength;
            
            musclelength = [muscle.muscle_length_profile(1); muscle.muscle_length_profile];
            
            
            
            dt = .54e-3;
            %dt = 2e-3;
            a = dt*ks/b;
            c = (1-a*(1+(kp/ks)));
            
            if abs(c) > 1
                error(['The passive tension solution for muscle ',muscle.muscle_name,' will oscillate to infinity.'])
            end
            
            % Solution based on notes from 3-5-20
            for i = 1:length(musclelength)
                if i == 1
                    delx(i,1) = musclelength(i)-Lr;
                    mV(i,1) = 0;
                    T(i,1) = ((ks*kp)/b)*delx(i,1)*dt;
                    Tdot(i,1) = T(i,1)/dt;
                elseif i == 2
                    delx(i,1) = musclelength(i)-Lr;
                    mV(i,1) = 0;
                    T(i,1) = c*T(i-1,1)+a*kp*delx(i,1)+a*b*mV(i,1);
                    Tdot(i,1) = (T(i,1)-T(i-1,1))/dt;
                else
                    delx(i,1) = musclelength(i)-Lr;
                    mV(i,1) = (delx(i,1)-delx(i-1,1))/dt;
                    T(i,1) = c*T(i-1,1)+a*kp*delx(i,1)+a*b*mV(i,1);
                    Tdot(i,1) = (T(i,1)-T(i-1,1))/dt;
                end
            end
            
            T(T<0) = 0;
            Tlen = min([length(muscleprofile) length(T)]);
            
             T_out(:,1) = muscleprofile(1:Tlen);
            T_out(:,2) = T(1:Tlen);
            
            %[~,dd] = findpeaks(-L_pass);
            tri = delaunayn(unique(T_out(:,1)));
            soughtindex = dsearchn(T_out(:,1),tri,length_find);
            T_find = T(soughtindex,1);
            
            return
            
            if (T(1)-T(end)) < 0
                    T_out(:,1) = L_pass;
                    T_out(:,2) = T;
            else

                avgblocks = 7;
                %converted to coefficients for the filtfilt function
                coeffblocks = ones(1,avgblocks)/avgblocks;

                Tsmooth = [T;flip(T);T;flip(T)];

                %For each joint, do the averaging of 10 data points many (originally set to
                %50) times

                for i=1:20
                    Tsmooth(:,1) = filtfilt(coeffblocks,1,Tsmooth(:,1));
                end
                Tsmooth = rescale(Tsmooth,min([T;flip(T);T;flip(T)]),max([T;flip(T);T;flip(T)]));
                Ts = Tsmooth(200:299,1);
                bb = T-Ts;
                cc = abs(gradient(bb));
                ccnorm = (cc-min(cc))/(max(cc)-min(cc));
                dd = gradient(ccnorm);
                [~,LOCS] = findpeaks(dd);
                gg = [Ts(1:LOCS(1));T(LOCS(1):end)];
                T_out(:,1) = L_pass;
                T_out(:,2) = gg(2:end);
                
%                 Tsmooth = rescale(Tsmooth,min([T;flip(T);T;flip(T)]),max([T;flip(T);T;flip(T)]));
%                 T_out(:,1) = L_pass;
%                 T_out(:,2) = Tsmooth(200:299,1);
            end
        end
        %% Function: lineInstersect3D: Find the closest point of intersection betweeen two vectors
        function P_intersect = lineIntersect3D(obj,PA,PB)
            % Find intersection point of lines in 3D space, in the least squares sense.
            % PA :          Nx3-matrix containing starting point of N lines
            % PB :          Nx3-matrix containing end point of N lines
            % P_Intersect : Best intersection point of the N lines, in least squares sense.
            % distances   : Distances from intersection point to the input lines
            % Anders Eikenes, 2012
            %from: https://www.mathworks.com/matlabcentral/fileexchange/37192-intersection-point-of-lines-in-3d-space

            Si = PB - PA; %N lines described as vectors
            ni = Si ./ (sqrt(sum(Si.^2,2))*ones(1,3)); %Normalize vectors
            nx = ni(:,1); ny = ni(:,2); nz = ni(:,3);
%             SXX = sum(nx.^2-1);
%             SYY = sum(ny.^2-1);
%             SZZ = sum(nz.^2-1);
%             SXY = sum(nx.*ny);
%             SXZ = sum(nx.*nz);
%             SYZ = sum(ny.*nz);
%             S = [SXX SXY SXZ;SXY SYY SYZ;SXZ SYZ SZZ];
%             CX  = sum(PA(:,1).*(nx.^2-1) + PA(:,2).*(nx.*ny)  + PA(:,3).*(nx.*nz));
%             CY  = sum(PA(:,1).*(nx.*ny)  + PA(:,2).*(ny.^2-1) + PA(:,3).*(ny.*nz));
%             CZ  = sum(PA(:,1).*(nx.*nz)  + PA(:,2).*(ny.*nz)  + PA(:,3).*(nz.^2-1));
%             C   = [CX;CY;CZ];
%             P_intersect = (S\C)';
            
            SXX_long = sum(reshape(nx.^2-1,2,length(nx)/2))';
            SYY_long = sum(reshape(ny.^2-1,2,length(ny)/2))';
            SZZ_long = sum(reshape(nz.^2-1,2,length(nz)/2))';
            SXY_long = sum(reshape(nx.*ny,2,length(nx)/2))';
            SXZ_long = sum(reshape(nx.*nz,2,length(nx)/2))';
            SYZ_long = sum(reshape(ny.*nz,2,length(ny)/2))';
            S_long = [reshape(SXX_long,1,1,length(SXX_long)),reshape(SXY_long,1,1,length(SXY_long)),reshape(SXZ_long,1,1,length(SXZ_long));...
                      reshape(SXY_long,1,1,length(SXY_long)),reshape(SYY_long,1,1,length(SYY_long)),reshape(SYZ_long,1,1,length(SYZ_long));...
                      reshape(SXZ_long,1,1,length(SXZ_long)),reshape(SYZ_long,1,1,length(SYZ_long)),reshape(SZZ_long,1,1,length(SZZ_long))];
            CX_long = sum(reshape(PA(:,1).*(nx.^2-1) + PA(:,2).*(nx.*ny)  + PA(:,3).*(nx.*nz),2,size(PA,1)/2))';
            CY_long = sum(reshape(PA(:,1).*(nx.*ny)  + PA(:,2).*(ny.^2-1) + PA(:,3).*(ny.*nz),2,size(PA,1)/2))';
            CZ_long = sum(reshape(PA(:,1).*(nx.*nz)  + PA(:,2).*(ny.*nz)  + PA(:,3).*(nz.^2-1),2,size(PA,1)/2))';
            C_long = squeeze([reshape(CX_long,1,1,length(CX_long));reshape(CY_long,1,1,length(CY_long));reshape(CZ_long,1,1,length(CZ_long))]);
            P_intersect = zeros(size(S_long,3),3);
            if all(S_long,'all')
                A = decomposition(S_long(:,:,1));
                P_intersect = (A\C_long)';
            else
                for ii = 1:size(S_long,3)
                    P_intersect(ii,:) = mldivide(S_long(:,:,ii),C_long(:,ii))';
                end
            end

%             if nargout>1
%                 N = size(PA,1);
%                 distances=zeros(N,1);
%                 for i=1:N %This is faster:
%                 ui=(P_intersect-PA(i,:))*Si(i,:)'/(Si(i,:)*Si(i,:)');
%                 distances(i)=norm(P_intersect-PA(i,:)-ui*Si(i,:));
%                 end
%             end
        end
        %% Function: Joint Axis Profile
        function [add_axis,R_axis] = joint_axis_profile(obj,jointMotion,joint)
            add_axis = zeros(3,size(jointMotion,1));
            R_axis = add_axis;
            for i =1:length(jointMotion)
                a = axis_angle_rotation(obj,jointMotion(i,1),obj.joint_obj{1}.uu_joint);
                b = axis_angle_rotation(obj,jointMotion(i,2),obj.joint_obj{2}.uu_joint);
                c = axis_angle_rotation(obj,jointMotion(i,3),obj.joint_obj{3}.uu_joint);
                
                if contains(joint.name,'Knee')
                    add_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR*[1;0;0])+joint.sim_position_profile(i,:)'*1000;
                    R_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR*[0;1;0])+joint.sim_position_profile(i,:)'*1000;
                elseif contains(joint.name,'Ankle')
                    add_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR*c*obj.body_obj{4}.CR*[1;0;0])+joint.sim_position_profile(i,:)'*1000;
                    R_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR*c*obj.body_obj{4}.CR*[0;1;0])+joint.sim_position_profile(i,:)'*1000;
                elseif contains(joint.name,'Hip')
                    add_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*[1;0;0])+joint.sim_position_profile(i,:)'*1000;
                    R_axis(:,i) = 10*(obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*[0;1;0])+joint.sim_position_profile(i,:)'*1000;
                else
                    keyboard
                end
                %joint_axis_w(:,i) = ankle_R_axis(:,i)*10+joint.sim_position_profile(i,:)'*1000;
            end
%             joint_axis_w(:,1:2000:20000);
%             jointMotion(1:2000:20000,1:3)';
        end
        %% Function: Find Muscle Index
        function number = find_muscle_index(obj,muscle)
            number = 0;
            if contains(muscle(end-2:end),'\n')
                if contains(muscle,'origin')
                    cutoff = strfind(muscle,'origin');
                else
                    if contains(muscle,'via')
                        cutoff = strfind(muscle,'via');
                    else
                        if contains(muscle,'ins')
                            cutoff = strfind(muscle,'ins');
                        else
                            keyboard
                        end
                    end
                end
                muscle = muscle(3:cutoff-1);    
            end
            for i = 1:length(obj.musc_obj)
                if contains(obj.musc_obj{i}.muscle_name,muscle)
                    number = i;
                    break
                end
            end
        end
        %% Function: Change Plot View with Input Angle
        function [az,el]=normalToAzimuthElevationDEG(obj,x,y,z,applyView)
            if nargin < 3
                applyView = 0;
            end
            if length(x)>1
                v         = x;
                if nargin > 1
                    applyView = y;
                end
                x=v(1);
                y=v(2);
                z=v(3);
            end
            if x==0 && y==0
                x =eps;
            end
            vNorm = sqrt(x^2+y^2+z^2);
            x=x/vNorm;
            y=y/vNorm;
            z=z/vNorm;
            az = (180/pi)*asin(x/sqrt(x^2+y^2));
            el = (180/pi)*asin(z);
            if applyView
                thereIsAnOpenFig = ~isempty(findall(0,'Type','Figure'));
                if thereIsAnOpenFig
                    axis equal
                    view([az,el]);
                    %the normal to the plane should now appear as a point
                    plot3([0,x],[0,y],[0,z],'linewidth',3)
                end
            end
        end
        %% Function: Compute the Spatial Manipulator Jacobian 
        function [Jac,p_foot] = compute_jac(obj,theta,thetaindex)
            %To clarify, this Jacobian is not the classical numerical jacobian
            %(differential of a function for each generalized coordinate). This
            %is the spatial manipulator Jacobian (see "A Mathematical Introduction to Robotic Manipulation" Murray, Li, and Sastry 1994)
            if length(theta) == length(obj.body_obj) - 1
                if size(theta,2) > size(theta,1)
                    theta = theta';
                end
                theta = [0;theta];
            elseif length(theta) == length(obj.body_obj)
                %do nothing
            elseif isempty(theta)
                theta = zeros(length(obj.body_obj),1);
            else
                disp('Orientation vector is not proper length. Should be a nx1 vector where n=num_bodies (first element is time)')
                p_foot = -1;
                Jac = -1;
                return
            end
            
            %thetaindex = find(obj.theta_motion(:,2)==theta(3),1);

            r_N = zeros(3,length(obj.body_obj));
            j_N = zeros(3,length(obj.body_obj)-1);
            
            for i = 2:length(obj.body_obj)
                %pos_bodies should be counted from 2 to end. Each column is that body's
                %position relative to the first body (which is 0s for the first).
                %Second body's position is given in the frame of the first.
                %r_N is the "World Position" in Animatlab of each body
                r_N(:,i) = obj.body_obj{i}.position_w_profile(thetaindex,:)';
                %Similarly, j_N is the "World Position" in Animatlab of each joint. 
                j_N(:,i-1) = obj.joint_obj{i-1,1}.sim_position_profile(thetaindex,:)';
            end
                
            %Axes of joint rotation
            omega = zeros(3,length(obj.body_obj)-1);

            %The jacobian matrix
            Jac = zeros(6,length(obj.body_obj)-1);

            %Use the local position of the toe.
            toe_pos = [20.733;-6.866;-1.867]/1000;

            a = axis_angle_rotation(obj,theta(2),obj.joint_obj{1}.uu_joint);
            b = axis_angle_rotation(obj,theta(3),obj.joint_obj{2}.uu_joint);
            c = axis_angle_rotation(obj,theta(4),obj.joint_obj{3}.uu_joint);

            foot_vec = r_N(:,4) + obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR*c*obj.body_obj{4}.CR*toe_pos;
           
            %This adds the starting world position of the toe (not after
            %theta rotation)
            j_N(:,end+1) = foot_vec;
            % Stores the world positions of the joint minus the hip joint
            p_st = j_N - repmat(j_N(:,1),1,size(j_N,2));
            %Hip joint
            leg_attach_pt = j_N(:,1);
            
            for i=1:length(obj.body_obj)-1
                omega(:,i) = obj.joint_obj{i}.uuw_joint(:,thetaindex);
                if strcmp(obj.joint_obj{i}.type,'Hinge')
                    %This is the joint twist. Crossing the axis of rotation
                    %with the member being rotated provides the twist
                    %Joint twists make up the top half of the spatial
                    %manipulator jacobian                    
                    Jac(1:3,i) = -cross(omega(:,i),j_N(:,i));
                    Jac(4:6,i) = omega(:,i);
                elseif strcmp(obj.joint_obj{i}.type,'Prismatic')
                    Jac(1:3,i) = omega(:,i);
                else
                    disp('Unidentified joint type.')
                end
            end
            
            p_foot = j_N(:,end)';
        end
        %% Function: Relevant Muscles
        function rmMat = find_relevant_muscles(obj)
            rmMat = cell(1,3);
            for i = 1:length(obj.musc_obj)
                attachment_bodies = cell2mat(obj.musc_obj{i}.pos_attachments(:,3));
                if attachment_bodies(1) == 1
                    rmMat{1} = [rmMat{1};i];
                end
                if attachment_bodies(end) == 3 || (attachment_bodies(1) == 2 && attachment_bodies(end) == 4)
                    rmMat{2} = [rmMat{2};i];
                end
                if attachment_bodies(end) == 4
                    rmMat{3} = [rmMat{3};i];
                end
            end
        end
        %% Function: Compute LOAD TORQUE
        function [lt_walk,load_torque] = compute_load_torques(obj,toplot)

                %Using a list of GRF at the foot, compute the maximum torques
                %that the joints should be able to apply.

                [beg,ennd,~] = find_step_indices(obj);
                [pks,locs] = findpeaks(obj.theta_motion(:,3));
                beg = locs(2); mid = locs(3); ennd = locs(4);
%                 beg = 654;
%                 ennd = 1676;

                %Midstance, ToeOff, MidSwing, ToeContact
                load([pwd,'\Data\MuirGRFData.mat'])
                m = size(obj.theta_motion(beg:ennd,:),1);
                n = length(VerticalNoSwing);
                animalmass = (obj.body_obj{1}.mass/44.4)*.30112;

                %Make kinematic swing data as long as stance data and put it all together
                lateralF = -interp1(1:n,LateralNoSwing,linspace(1,n,m))*animalmass;
                verticalF = interp1(1:n,VerticalNoSwing,linspace(1,n,m))*animalmass;
                propulsiveF = interp1(1:n,PropulsiveNoSwing,linspace(1,n,m))*animalmass;

                    forcesmooth = [propulsiveF,propulsiveF,propulsiveF;verticalF,verticalF,verticalF;lateralF,lateralF,lateralF]';
                    forcesmooth = smoothdata(forcesmooth,'gaussian',100);

                    %sizer = floor(size(forcesmooth,1)/6);
                    temp = find(diff(forcesmooth(:,1)==0)==-1);
                    forcesmooth = forcesmooth(temp(1):temp(2),:);
                    
                    %This needs to be rotated to match the Animatlab xyz environment. This is relatively simple at the moment, we just need to rotate planar
                    %forces by about 16 degrees clockwise along the global Z axis.
                    theta = -16.3*(pi/180);
                    rotmat = [cos(theta) -sin(theta) 0;...
                              sin(theta) cos(theta) 0;...
                              0 0 1];
                    forcesmooth = rotmat*forcesmooth';
                    
                    load_torque = zeros(size(lateralF,2),3);
                    foot_position = load_torque;
                    
                    gWH = [[obj.body_obj{1}.CR [0;0;0]];[zeros(1,3) 1]];
                    gHF = [[obj.body_obj{2}.CR obj.body_obj{2}.position'];[zeros(1,3) 1]];
                    gFT = [[obj.body_obj{3}.CR obj.body_obj{3}.position'];[zeros(1,3) 1]];
                    gTX = [[obj.body_obj{4}.CR obj.body_obj{4}.position'];[zeros(1,3) 1]];
                for i=beg:ennd
                    pointer = i-beg+1;

                    current_config = obj.theta_motion(i,:)';

                    %Input the joint angles that we just computed, and
                    %insert them into the configuration of the leg.
                    [J,foot_position(pointer,:)] = obj.compute_jac(current_config,i);    

                    %Flip signs because this is the force the toe is pushing against the ground. This analysis is effector focused so motion of the
                    %effector, forces made by the effector, etc.
%                     force_vec = [lateralF(pointer);verticalF(pointer);propulsiveF(pointer)];
                    force_vec = forcesmooth(:,pointer);

                    %Calculate J'*F to find the joint torques. 
                    %For reference, consult Sastry 1994 p. 121 (Eq 3.60)
                    %foot_wrench_body = [force_vec(:,k);cross(foot_position(:,i),force_vec(:,k))];
                    foot_wrench = [force_vec;0;0;0];
                    
%                     a = [[axis_angle_rotation(obj,current_config(1),obj.joint_obj{1}.uu_joint);zeros(1,3)],[zeros(3,1);1]];
%                     b = [[axis_angle_rotation(obj,current_config(2),obj.joint_obj{2}.uu_joint);zeros(1,3)],[zeros(3,1);1]];
%                     c = [[axis_angle_rotation(obj,current_config(3),obj.joint_obj{3}.uu_joint);zeros(1,3)],[zeros(3,1);1]];
                    a = [[axis_angle_rotation(obj,current_config(1),obj.joint_obj{1}.uuw_joint(:,i));zeros(1,3)],[zeros(3,1);1]];
                    b = [[axis_angle_rotation(obj,current_config(2),obj.joint_obj{2}.uuw_joint(:,i));zeros(1,3)],[zeros(3,1);1]];
                    c = [[axis_angle_rotation(obj,current_config(3),obj.joint_obj{3}.uuw_joint(:,i));zeros(1,3)],[zeros(3,1);1]];
                    gTot = gWH*a*gHF*b*gFT*c*gTX;
                    gTotR = obj.body_obj{1}.CR*axis_angle_rotation(obj,current_config(1),obj.joint_obj{1}.uuw_joint(:,i))*...
                        obj.body_obj{2}.CR*axis_angle_rotation(obj,current_config(2),obj.joint_obj{2}.uuw_joint(:,i))*...
                        obj.body_obj{3}.CR*axis_angle_rotation(obj,current_config(3),obj.joint_obj{3}.uuw_joint(:,i))*obj.body_obj{4}.CR;
                    gTotP = obj.musc_obj{20}.pos_attachments{5,4}(i,:)';
                    
%                     gTotR = gTot(1:3,1:3);
%                     gTotP = gTot(1:3,4);
                    
                    invAd = [[gTotR';zeros(3,3)],[-gTotR'*wedge(gTotP);gTotR']];
                    Js = invAd*J;
                    
                    current_torques = Js'*foot_wrench;
                    
                    %The torques in each joint as generated only by the force of the ground on the foot
                    load_torque(pointer,:) = current_torques';
                end
                    load_torque = interp1(1:length(load_torque),load_torque,linspace(1,length(load_torque),length(obj.theta_motion)));
                % Fit the load torque to the walking waveform. Find a way to do this automatically. In the meantime, just say explcitly when stance
                % and swing happen in the walking waveform (obj.theta_motion)
                % translocs represent the beginning of motion, the peaks and troughs of the hip motion, and the end of motion
                translocs = [1,locs(2:end)',length(obj.theta_motion)];
                shifter = floor((translocs(3)-translocs(2))*.95);
                translocs = translocs - [0,0,shifter,0,shifter,0,0];
                translocs = [1 657 853 1682 1878 2706 length(obj.theta_motion)];
                translocs = [1 584 665 1611 1694 2638 2709 length(obj.theta_motion)];
                % Find translocs
                [~,hlocs] = findpeaks(obj.theta_motion(:,1));
                [~,alocs] = findpeaks(obj.theta_motion(:,3));
                translocs = [1 hlocs(2) alocs(2) hlocs(3) alocs(4) hlocs(4) alocs(5) length(obj.theta_motion)];
                loadEnd = find(diff(load_torque(:,3)==0)==1,1,'first');
                stance = load_torque(1:loadEnd-1,:); swing = load_torque(loadEnd:end,:);
                lt_walk = [];
                for ii = 3:2:length(translocs)
                    st = interp1(1:length(stance),stance,linspace(1,length(stance),translocs(ii-1)-translocs(ii-2)));
                    sw = interp1(1:length(swing),swing,linspace(1,length(swing),translocs(ii)-translocs(ii-1)));
                    lt_walk = [lt_walk;st;sw];
                end
                %st = interp1(1:length(stance),stance,linspace(1,length(stance),length(obj.theta_motion)-length(lt_walk)));
                temp  =length(obj.theta_motion)-length(lt_walk);
                lt_walk = [lt_walk;st(1:temp,:)];
                
                %Plotting
                if toplot
                    figure('Position',[962,2,958,994]);
                        ylimms1(1) = 1.1*min(min(forcesmooth));
                        ylimms1(2) = 1.1*max(max(forcesmooth));
                        ylimms2(1) = 1.1*min(min(1000*load_torque));
                        ylimms2(2) = 1.1*max(max(1000*load_torque));
                        subA = subplot(3,1,1);
                            plot(100*linspace(0,1,size(obj.theta_motion(locs(3):locs(5),:),1)),obj.theta_motion(locs(3):locs(5),:).*(180/pi)+[98.4373,102.226,116.2473],'LineWidth',2)
                            title('Joint Angles for Stance and Swing')
                            ylabel('Angle (deg)')
                            xlabel('% Stride')
                            legend({'Hip','Knee','Ankle'},'Location','eastoutside')
                            grid on
                        subB = subplot(3,1,2);
                            plot(100*linspace(0,1,size(forcesmooth',1)),forcesmooth','LineWidth',2)
                            %hold on
                            %plot(linspace(0,1,size(bb,1)),bb)
                            title('Ground Reaction Forces')
                            xlabel('% Stride')
                            ylabel('Force (N)')
                            legend({'Propulsive','Vertical','Lateral'},'Location','eastoutside')
                            ylim(ylimms1)
                            grid on
                        subC = subplot(3,1,3);
                            plot(100*linspace(0,1,size(load_torque,1)),1000*load_torque','LineWidth',2)
                            title('Load Torques')
                            xlabel('% Stride')
                            ylabel('Torque (mN-m)')
                            legend({'Hip','Knee','Ankle'},'Location','eastoutside')
                            ylim(ylimms2)
                            grid on
                            
                            pause(.5)
                            subApos = get(subA,'Position');
                            subBpos = get(subB,'Position');
                            set(subB,'Position',[subApos(1) subBpos(2) subApos(3) subApos(4)]);

                            ax = findobj(gcf,'Type','Axes');
                            for i=1:length(ax)
                                set(ax(i).XLabel,'FontSize',14)
                                set(ax(i).YLabel,'FontSize',14)
                                set(ax(i).Title,'FontSize',16)
                            end
                         
                        %saveas(gcf,[pwd,'\OutputFigures\Images\computer_active_joint_torque\','total_load_torque_',datestr(datetime('now'),'yyyymmdd'),'.png'])
                end
        end
        %% Function: Compute BODY TORQUE
        function grav_torque = compute_grav_torques(obj,toplot)
            [beg,ennd,~] = find_step_indices(obj);
            %Higher div, the longer the program takes.
            div = 500;
            xx = obj.sampling_vector;
            
            gravDir = [0,-9.8,0]; % Gravity is along the Y-axis at -9.8 m/s^2
            
            bodyMass = zeros(1,length(obj.body_obj));
            for ii = 1:length(obj.body_obj)
                % Body mass in kg
                bodyMass(ii) = obj.body_obj{ii}.mass./1000;
            end
                
            grav_torque = zeros(length(xx),length(obj.joint_obj));
            for jj = 1:length(obj.joint_obj)
                [bodyArm] = compute_gravity_moment_arm(obj,obj.body_obj{jj+1},jj,0);
                bodyTorque = (bodyArm(:,3)./1000).*bodyMass(jj+1).*9.8;
                if jj==3
                    nJointArm = zeros(size(bodyArm));
                else
                    [nJointArm] = compute_gravity_moment_arm(obj,obj.joint_obj{jj+1},jj,0);
                end
                nJointTorque = (nJointArm(:,3)./1000).*sum(bodyMass(jj+2:end)).*9.8;
                grav_torque(:,jj) = bodyTorque'+nJointTorque';
            end
            grav_torque = -grav_torque;
            
            clearTheta = obj.theta_motion.*(180/pi)+[98.4373 102.226 116.2473];
            if toplot
                figure('name','GravTorques','Position',[962,2,958,994]);
                subplot(2,1,1)
                    plot(obj.theta_motion_time(xx),-grav_torque,'LineWidth',3)
                    legend({'Hip';'Knee';'Ankle'})
                    ylabel('Joint Torque (mNm)')
                    xlabel('Time (s)')
                    xlim([0,max(obj.theta_motion_time(xx))])
                    title('Gravitational Torques from Body Segments','FontSize',16)
                subplot(2,1,2)
                    plot(obj.theta_motion_time(xx),clearTheta,'LineWidth',3)
                    legend({'Hip';'Knee';'Ankle'})
                    ylabel('Joint Motion (XX)')
                    xlabel('Time (s)')
                    xlim([0,max(obj.theta_motion_time(xx))])
                    title('Joint Motion','FontSize',16)
            end
        end
        %% Function: Compute MUSCLE TORQUE: Compute the passive joint torque from muscles for a joint over an entire walking cycle
        function [passive_joint_torque,passive_joint_motion,passive_muscle_torque] = compute_passive_muscle_torque(obj)
            num_muscles = length(obj.musc_obj);
            passive_joint_motion = obj.theta_motion.*(180/pi)+[98.4373,102.226,116.2473];
            passive_joint_torque = zeros(length(obj.theta_motion),3);
            relevant_muscles = find_relevant_muscles(obj);   
            passive_muscle_torque = zeros(num_muscles,length(obj.theta_motion),3);

            for k = 1:3   
                [moment_output] = compute_joint_moment_arms(obj,k,1);
                Tpass_profile = zeros(size(moment_output));
                for i = 1:length(relevant_muscles{k})
                    muscle_num = relevant_muscles{k}(i);
                    musc = obj.musc_obj{muscle_num};
                    % This will trim the rapid changes in passive tension at the beginning and end of the sim, which are artifacts
                    % of moving the leg into and out of position. Doing this avoids huge values for the passive torque, which can affect
                    % the force optimization process. Feel free to remove this line and things will still work (set pt = musc.passive_tension)
                    lenDat = length(musc.passive_tension);
                    pTemp = [repmat(musc.passive_tension(34),33,1);musc.passive_tension(34:lenDat-9);repmat(musc.passive_tension(lenDat-9),9,1)];
                    pt = smoothdata(pTemp,'movmedian',25);
                    Tpass_profile(muscle_num,:) = pt;
                    passive_muscle_torque(muscle_num,:,k) = pt'.*moment_output(muscle_num,:)./1000;
                end
                passive_joint_torque(:,k) = smoothdata(sum(passive_muscle_torque(:,:,k)),'loess',25);
            end
            plotter = 0;
            if plotter
                figure('Position',[962,2,958,994]);
                plot(obj.theta_motion_time,passive_joint_torque')
                title('Passive Joint Torque')
                xlabel('Time (s)')
                ylabel('Torque (Nm)')
                legend({'Hip','Knee','Ankle'})
                for k = 1:3
                    figure('Position',[962,2,958,994]);
                    count = 1;
                    CM = hsv(size(relevant_muscles{k},1));
                    legendcell = {};
                    for i = 1:length(relevant_muscles{k})
                        muscle_num = relevant_muscles{k}(i);
                        plot(obj.theta_motion_time(100:end-100),1000*passive_muscle_torque(muscle_num,100:end-100,k),'color',CM(count,:),'LineWidth',1.5)
                        legendcell{end+1} = obj.musc_obj{muscle_num}.muscle_name(4:end);
                        count = count + 1;
                        %title(obj.musc_obj{i}.muscle_name,'Interpreter','none')
                        hold on
                    end
                    title(['Passive Muscle Torque of the ',obj.joint_obj{k}.name(4:end)])
                    xlabel('Time (s)')
                    ylabel('Torque (mN-m)')
                    legend(legendcell,'Interpreter','none','Location','eastoutside')
                    %set(gcf,'Position',[500 500 900 500])
                    %saveas(gcf,['G:\My Drive\Rat\Optimizer\OutputFigures\passive_torque','_',obj.joint_obj{k}.name(4:end),'_',datestr(datetime('now'),'yyyymmdd'),'.png'])
                end
            end
        end
        %% Function: Compute INERTIAL TORQUE
        function [inertial_torque] = compute_inertial_torque(obj)
            % Find angular acceleration
            g = 9.81;
            lenDat = length(obj.theta_motion);
            jm = [repmat(obj.theta_motion(21,:),20,1);obj.theta_motion(21:lenDat-9,:);repmat(obj.theta_motion(lenDat-9,:),9,1)];
            jm = smoothdata(jm,'movmean',50);
            gjm = zeros(length(jm),3); gjm2 = gjm;
            for ii = 1:3
                gjm(:,ii) = smoothdata(gradient(jm(:,ii),obj.dt_motion),'movmean',100);
                gjm2(:,ii) = smoothdata(gradient(gjm(:,ii),obj.dt_motion),'movmean',100);
            end
            gjm2 = [repmat(gjm2(56,:),56,1);gjm2(56:lenDat-33,:);repmat(gjm2(lenDat-33,:),35,1)];
            moi2 = zeros(3,1);
            % Computing the moments of inertia
            % These are just approximations, assuming femur and tibia are solid rods and the foot is a box
            % Measurements for dimensions approximated by inspection using Animatlab
            % Notes on 6-2-20
%             MOI = zeros(3,3,3);
                m1 = obj.body_obj{2}.mass/1000;
                m2 = obj.body_obj{3}.mass/1000;
                m3 = obj.body_obj{4}.mass/1000;
%             % Femur 
%             L1 = obj.body_obj{2}.length; r1 = 0.049654*L1; MOI(1,1,1) = (1/12)*m1*L1^2; MOI(2,2,1) = (1/12)*m1*L1^2; MOI(3,3,1) = (1/2)*m1*r1^2;
%                 comlen1 = norm(com_pos_on_demand(obj,obj.theta_motion(1,:),2)' - obj.joint_obj{1}.sim_position_profile(1,:));
%                 moi2(1) = (1/3)*m1*L1^2+(m2+m3)*L1^2;
%             % Tibia
%             L2 = obj.body_obj{3}.length; r2 = 0.03061*L2; MOI(1,1,2) = sqrt(((1/12)*m2*L2^2)^2+((1/2)*m2*r2^2)^2); MOI(2,2,2) = sqrt(((1/12)*m2*L2^2)^2+((1/2)*m2*r2^2)^2); MOI(3,3,2) = (1/12)*m2*L2^2;
%                 comlen2 = norm(com_pos_on_demand(obj,obj.theta_motion(1,:),3)' - obj.joint_obj{2}.sim_position_profile(1,:));
%                 moi2(2) = (1/3)*m2*L2^2+m3*L1^2;
%             % Foot
%             d = obj.body_obj{4}.length; h = 0.13473*d; w = 0.44908*d; MOI(1,1,3) = (1/12)*m3*(d^2+h^2); MOI(2,2,3) = (1/12)*m3*(d^2+w^2); MOI(3,3,3) = (1/12)*m3*(h^2+w^2);
%                 comlen3 = norm(com_pos_on_demand(obj,obj.theta_motion(1,:),4)' - obj.joint_obj{3}.sim_position_profile(1,:));
%                 moi2(3) = (1/12)*m3*(h^2 + d^2)+m3*comlen3^2;
%             inertial_torque = -gjm2.*moi2';
            
            for ii = 1:length(obj.theta_motion)
                M = compute_mass_matrix(obj,jm(ii,:));
                %[C,M] = compute_coriolis_matrix(obj,jm(ii,:),gjm(ii,:));
%                     hippos = obj.joint_obj{1}.sim_position_profile(ii,2);
%                     fempos = com_pos_on_demand(obj,jm(ii,:),2)';
%                     h(ii,1) = fempos(2)-hippos;
%                     tibpos = com_pos_on_demand(obj,jm(ii,:),3)';
%                     h(ii,2) = tibpos(2)-hippos;
%                     footpos = com_pos_on_demand(obj,jm(ii,:),4)';
%                     h(ii,3) = footpos(2)-hippos;
%                     V(ii,:) = [m1*g,m2*g,m3*g].*h(ii,:);
                inertial_torque(ii,:) = -19.42*(M*gjm2(ii,:)')';
            end
        end
        %% Function: Compute ACTIVE TORQUE: Passive Joint Torque
        function  [total_joint_torques,passive_muscle_torque,inertial_torque] = compute_active_joint_torque(obj,toplot)
                        
            % Muscle torque
            [passive_muscle_torque,passive_joint_motion] = compute_passive_muscle_torque(obj);
            
            % Load Torque
            load_torque = compute_load_torques(obj,0);
            
            % Body Torque
            passive_grav_torque = compute_grav_torques(obj,0);
            
            % Inertial Torque
            inertial_torque = compute_inertial_torque(obj);
            
%             %
%             passive_muscle_torque = -passive_muscle_torque;
%             load_torque = -load_torque;
%             passive_grav_torque = -passive_grav_torque;
%             inertial_torque = -inertial_torque;
            
            load_torque_on = 1;
            if load_torque_on
                total_joint_torques = passive_grav_torque+passive_muscle_torque+inertial_torque+load_torque;

                %Plotting
                if toplot
                    scale = 1.2;
                    aa = figure('Position',[962,2,958,994]);
                    subplot(6,1,1)
                        plot(obj.theta_motion_time,passive_joint_motion','LineWidth',2)
                        title('Joint Motion'); ylabel('Angle (deg)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms1 = [min(passive_joint_motion,[],'all') max(passive_joint_motion,[],'all')];
                        ylim(ylimms1); xlim([0 max(obj.theta_motion_time)])
                    subplot(6,1,2)
                        plot(obj.theta_motion_time,1000*passive_muscle_torque,'LineWidth',2)
                        title('Passive Torque (Muscle)'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms2(1) = 1000*scale*min(passive_muscle_torque(50:end,:),[],'all');
                        ylimms2(2) = 1000*scale*max(passive_muscle_torque(50:end,:),[],'all');
                        ylim(ylimms2); xlim([0 max(obj.theta_motion_time)])
                    subplot(6,1,3)
                        plot(obj.theta_motion_time,1000.*passive_grav_torque,'LineWidth',2)
                        title('Body Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        xlim([0 max(obj.theta_motion_time)])
                    subplot(6,1,4)
                        plot(obj.theta_motion_time,1000.*inertial_torque,'LineWidth',2)
                        title('Inertial Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        xlim([0 max(obj.theta_motion_time)])
                    subplot(6,1,5)
                        plot(obj.theta_motion_time,1000.*load_torque,'LineWidth',2)
                        title('Load Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        xlim([0 max(obj.theta_motion_time)])
                    subplot(6,1,6)
                        plot(obj.theta_motion_time,1000.*total_joint_torques,'LineWidth',2)
                        title('Total Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms2(1) = 1000*scale*min(total_joint_torques(50:end,:),[],'all');
                        ylimms2(2) = 1000*scale*max(total_joint_torques(50:end,:),[],'all');
                        ylim(ylimms2); xlim([0 max(obj.theta_motion_time)])
                end
            else
                total_joint_torques = passive_grav_torque+passive_muscle_torque+inertial_torque;

                %Plotting
                if toplot
                    scale = 1.2;
                    aa = figure('Position',[962,2,958,994]);
                    subplot(5,1,1)
                        plot(obj.theta_motion_time,passive_joint_motion','LineWidth',2)
                        title('Joint Motion'); ylabel('Angle (deg)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms1 = [min(passive_joint_motion,[],'all') max(passive_joint_motion,[],'all')];
                        ylim(ylimms1); xlim([0 max(obj.theta_motion_time)])
                    subplot(5,1,2)
                        plot(obj.theta_motion_time,1000*passive_muscle_torque,'LineWidth',2)
                        title('Passive Torque (Muscle)'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms2(1) = 1000*scale*min(passive_muscle_torque(50:end,:),[],'all');
                        ylimms2(2) = 1000*scale*max(passive_muscle_torque(50:end,:),[],'all');
                        ylim(ylimms2); xlim([0 max(obj.theta_motion_time)])
                    subplot(5,1,3)
                        plot(obj.theta_motion_time,1000.*passive_grav_torque,'LineWidth',2)
                        title('Body Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        xlim([0 max(obj.theta_motion_time)])
                    subplot(5,1,4)
                        plot(obj.theta_motion_time,1000.*inertial_torque,'LineWidth',2)
                        title('Inertial Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        xlim([0 max(obj.theta_motion_time)])
                    subplot(5,1,5)
                        plot(obj.theta_motion_time,1000.*total_joint_torques,'LineWidth',2)
                        title('Total Torques'); ylabel('Torque (mN-m)'); xlabel('Time (s)')
                        legend('Hip','Knee','Ankle','Location','eastoutside')
                        ylimms2(1) = 1000*scale*min(total_joint_torques(50:end,:),[],'all');
                        ylimms2(2) = 1000*scale*max(total_joint_torques(50:end,:),[],'all');
                        ylim(ylimms2); xlim([0 max(obj.theta_motion_time)])
                end    
            end
        end
        %% Function: Plot Joint Axes Over Time
        function joint_axes_vis(obj)
            %This function is meant to help visualize each joint's uu_joint, uu_joint2, and uu_joint3. Use this to help make Ext/Flx,Abd/Add,Exr/InR motion
            %make sense. Remember: uu_joints are all defined in the distal body's coordinate system (knee->tibia, ankle->foot) 
            scale = 1000;
            i = 1;
            pausebutton = uicontrol('Style', 'ToggleButton', ...
                'Units',    'pixels', ...
                'Position', [125 400 60 20], ...
                'String',   'Pause', ...
                'Value',    0);
            stopbutton = uicontrol('Style', 'ToggleButton', ...
                'Units',    'pixels', ...
                'Position', [125 370 60 20], ...
                'String',   'Debug', ...
                'Value',    0);
            while i <= length(obj.theta_motion)
                %for i=1:10:length(obj.theta_motion)
%                 a = obj.joint_obj{1,1}.joint_rotmat_profile(:,:,i);
%                 b = obj.joint_obj{2,1}.joint_rotmat_profile(:,:,i);
%                 c = obj.joint_obj{3,1}.joint_rotmat_profile(:,:,i);
                femur = [scale*obj.joint_obj{1}.sim_position_profile(i,1),scale*obj.joint_obj{1}.sim_position_profile(i,2),scale*obj.joint_obj{1}.sim_position_profile(i,3);...
                scale*obj.joint_obj{2}.sim_position_profile(i,1),scale*obj.joint_obj{2}.sim_position_profile(i,2),scale*obj.joint_obj{2}.sim_position_profile(i,3)];
                tibia = [scale*obj.joint_obj{2}.sim_position_profile(i,1),scale*obj.joint_obj{2}.sim_position_profile(i,2),scale*obj.joint_obj{2}.sim_position_profile(i,3);...
                scale*obj.joint_obj{3}.sim_position_profile(i,1),scale*obj.joint_obj{3}.sim_position_profile(i,2),scale*obj.joint_obj{3}.sim_position_profile(i,3)];
                foot = scale*[obj.joint_obj{3}.sim_position_profile(i,1),obj.joint_obj{3}.sim_position_profile(i,2),obj.joint_obj{3}.sim_position_profile(i,3);...
                            obj.musc_obj{20, 1}.pos_attachments{5,4}(i,:)];
                hippos = scale*obj.joint_obj{1}.sim_position_profile(i,:);
                kneepos = scale*obj.joint_obj{2}.sim_position_profile(i,:);
                anklepos = scale*obj.joint_obj{3}.sim_position_profile(i,:);
                %Plot the hip
                    plot3(scale*obj.joint_obj{1}.sim_position_profile(i,1),scale*obj.joint_obj{1}.sim_position_profile(i,2),scale*obj.joint_obj{1}.sim_position_profile(i,3),'kp')
                    hold on
                    %plot hip uu_joint
                        hipuu  = [hippos;hippos+(obj.joint_obj{1}.uuw_joint)'*10];
                        hipuu2 = [hippos;hippos+(obj.joint_obj{1}.uuw_joint2(:,i))'*10];
                        hipuu3 = [hippos;hippos+(obj.joint_obj{1}.uuw_joint3(:,i))'*10];
%                         hipuu  = [hippos;hippos+(obj.CR_bodies(:,:,1)*obj.joint_obj{1}.uu_joint)'*10];
%                         hipuu2 = [hippos;hippos+(obj.CR_bodies(:,:,1)*obj.joint_obj{1}.uu_joint2(:,i))'*10];
%                         hipuu3 = [hippos;hippos+(obj.CR_bodies(:,:,1)*obj.joint_obj{1}.uu_joint3(:,i))'*10];
                        plot3(hipuu(:,1),hipuu(:,2),hipuu(:,3),'r')
                        %plot hip uu_joint2
                        plot3(hipuu2(:,1),hipuu2(:,2),hipuu2(:,3),'g')
                        %plot hip uu_joint3
                        plot3(hipuu3(:,1),hipuu3(:,2),hipuu3(:,3),'b')
                %Plot the knee
                    plot3(scale*obj.joint_obj{2}.sim_position_profile(i,1),scale*obj.joint_obj{2}.sim_position_profile(i,2),scale*obj.joint_obj{2}.sim_position_profile(i,3),'ks')
                    hold on
                    %plot knee uu_joint
                        kneeuu  = [kneepos;kneepos+(obj.joint_obj{2}.uuw_joint)'*10];
                        kneeuu2 = [kneepos;kneepos+(obj.joint_obj{2}.uuw_joint2(:,i))'*10];
                        kneeuu3 = [kneepos;kneepos+(obj.joint_obj{2}.uuw_joint3(:,i))'*10];
%                         kneeuu  = [kneepos;kneepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.joint_obj{2}.uu_joint)'*10];
%                         kneeuu2 = [kneepos;kneepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.joint_obj{2}.uu_joint2(:,i))'*10];
%                         kneeuu3 = [kneepos;kneepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.joint_obj{2}.uu_joint3(:,i))'*10];
                        plot3(kneeuu(:,1),kneeuu(:,2),kneeuu(:,3),'r')
                        %plot knee uu_joint2
                        plot3(kneeuu2(:,1),kneeuu2(:,2),kneeuu2(:,3),'g')
                        %plot knee uu_joint3
                        plot3(kneeuu3(:,1),kneeuu3(:,2),kneeuu3(:,3),'b')  
                %Plot the femur
                plot3(femur(:,1),femur(:,2),femur(:,3),'k--')
                %Plot the tibia
                plot3(tibia(:,1),tibia(:,2),tibia(:,3),'k-.')
                %Plot the ankle
                    plot3(scale*obj.joint_obj{3}.sim_position_profile(i,1),scale*obj.joint_obj{3}.sim_position_profile(i,2),scale*obj.joint_obj{3}.sim_position_profile(i,3),'kd')
                    hold on
                    %plot ankle uu_joint
                        ankleuu  = [anklepos;anklepos+(obj.joint_obj{3}.uuw_joint)'*10];
                        ankleuu2 = [anklepos;anklepos+(obj.joint_obj{3}.uuw_joint2(:,i))'*10];
                        ankleuu3 = [anklepos;anklepos+(obj.joint_obj{3}.uuw_joint3(:,i))'*10];
%                         ankleuu  = [anklepos;anklepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.joint_obj{3}.uu_joint)'*10];
%                         ankleuu2 = [anklepos;anklepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.joint_obj{3}.uu_joint2(:,i))'*10];
%                         ankleuu3 = [anklepos;anklepos+(obj.CR_bodies(:,:,1)*obj.CR_bodies(:,:,2)*obj.CR_bodies(:,:,3)*obj.joint_obj{3}.uu_joint3(:,i))'*10];
                        plot3(ankleuu(:,1),ankleuu(:,2),ankleuu(:,3),'r')
                        %plot ankle uu_joint2
                        plot3(ankleuu2(:,1),ankleuu2(:,2),ankleuu2(:,3),'g')
                        %plot ankle uu_joint3
                        plot3(ankleuu3(:,1),ankleuu3(:,2),ankleuu3(:,3),'b')  
                %Plot the foot
                plot3(foot(:,1),foot(:,2),foot(:,3),'k-.')
                grid on
                %Oblique View
%                 view([-32 31])
                view([0 90])
                xlim([-50 50])
                ylim([-80 10])
                zlim([0 30])
                axeslengths = [max(xlim)-min(xlim);max(ylim)-min(ylim);max(zlim)-min(zlim)];
                normedaxes = axeslengths/norm(axeslengths);
                pbaspect(normedaxes)
                set(gcf,'Position',[500 100 1300 800]);
                figure(1)
                hold off
                i = i + 10;
                pause_value = get(pausebutton, 'Value');
                stop_value = get(stopbutton, 'Value');
                pause(.001)
                if stop_value
                    keyboard
                end
                while pause_value
%                     pause
%                     set(buttonH,'Value',0)
%                     button_value = 0;
                    pause(.001)
                    pause_value = get(pausebutton, 'Value');
                    if get(stopbutton, 'Value')
                        keyboard
                    end
                end  
            end
            %end
        end
        %% Function: Plot and Classify Muscles by Muscle Lengths vs. Joint Angles
        function plot_muscles_by_joint(obj,jointnum,fullstep)
            relevant_muscles = [];
            num_muscles = length(obj.musc_obj);
            
            %classifyby = 1 is by muscle length
            %classifyby = 2 is by momentarm
            classifyby = 2;
            
            if fullstep ~= 0 && fullstep ~=1
                fullstep = input('Please specify whether you want to plot a full step of motion (1) or want to see the joint range of motion (0)');
                if fullstep ~= 0 && fullstep ~=1
                    fullstep = 1;
                end
            end
            
            for i = 1:num_muscles
                attachment_bodies = cell2mat(obj.musc_obj{i}.pos_attachments(:,3));
                if jointnum == 1 
                    if attachment_bodies(1) == 1
                        relevant_muscles = [relevant_muscles;i];
                    end
                elseif jointnum == 2
                    if attachment_bodies(end) == 3
                        relevant_muscles = [relevant_muscles;i];
                    end
                elseif jointnum == 3
                    if attachment_bodies(end) == 4
                        relevant_muscles = [relevant_muscles;i];
                    end
                end
            end
            
            if fullstep
                [beg,ennd,~] = find_step_indices(obj);
            else
                jointprofile = obj.joint_obj{jointnum}.rec_angle_profile;
                [~,maxlocs] = findpeaks(jointprofile);
                [~,minlocs] = findpeaks(-jointprofile);
                a3 = unique([maxlocs;minlocs]);
                a4 = jointprofile(a3);
                for i = 1:length(a4)-1
                    range(i,1) = abs(a4(i+1,1)-a4(i,1));
                end
                [~,dd] = max(range);
                dd = dd(1);
                beg = a3(dd);
                ennd = a3(dd+1);
            end

            fig1count = 0;
            fig2count = 0;
            fig3count = 0;
            
            %Change subplots to 4 to add an extra plot that shows the signs of the muscle derivative. Used to classify the muscles.
            subplots = 4;
            figure
            for i = 1:size(relevant_muscles,1)
                if classifyby == 2
                    [moment_arm_profile,mom_arm_inst_length] = compute_muscle_moment_arm(obj,obj.musc_obj{relevant_muscles(i)},1,jointnum,0,0);
                    muscle_length = moment_arm_profile(:,3);
                    time = moment_arm_profile(:,1);
                    mnorm = muscle_length/(max(muscle_length) - min(muscle_length));
                else
                    muscle_length = obj.musc_obj{relevant_muscles(i)}.muscle_length_profile(beg:ennd);
                    time = obj.joint_obj{jointnum}.rec_angle_time(beg:ennd);
                    mnorm = (muscle_length - min(muscle_length)) / (max(muscle_length) - min(muscle_length));
                end
                jointangle = obj.joint_obj{jointnum}.rec_angle_profile(beg:ennd)*(180/pi);
                tnorm = (time - min(time)) / (max(time) - min(time));
                if fullstep
                    dy=smooth(diff(mnorm)./diff(tnorm));
                else
                    dy=smooth(diff(mnorm)./diff(jointangle));
                end
                a = size(dy);
                dy = dy(~any(isnan(dy) | isinf(dy) | dy==0,2),:);
                b = size(dy);
                dy = sign(dy);
                jointangle = jointangle(2+(a-b):end);
                mnorm = mnorm(2+(a-b):end);
                tnorm = tnorm(2+(a-b):end);         
                dd = diff(dy);
                dd = dd(dd~=0);
                transitions = find(diff(dy));
                switches = size(dd,1);
                if size(unique(dy),1) > 1 
                    [c,~]=histogram(dy,unique(dy));
                end
                
                if subplots == 4
                    subplot(subplots,1,4)
                    if fullstep
                        plot(tnorm,dy)
                        title(['Slope switches =',num2str(switches),' firsttrans =',num2str(dd(1))]);
                    else
                        plot(jointangle,dy)
                        title([num2str(imp),' numpos=',num2str(bb),' numneg=',num2str(cc)]);
                    end
                end
                if classifyby == 2
                    if fullstep
                        if jointnum == 1
                            if dd(1) == -2
                                figurecount = 1;
                                fig1count = fig1count + 1;
                                legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                if mean(dy(1:floor(size(dy,1)/2))) > 0
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                else
                                    figurecount = 3;
                                    fig3count = fig3count + 1;
                                    legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                end
                            end
%                             if switches > 1 && abs(mean(dy(1:floor(size(dy,1)/2)))) < .9
%                                 figurecount = 3;
%                                 fig3count = fig3count + 1;
%                                 legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
%                             else
%                                 if mean(dy(1:floor(size(dy,1)/2))) > 0
%                                     figurecount = 2;
%                                     fig2count = fig2count + 1;
%                                     legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
%                                 else
%                                     figurecount = 1;
%                                     fig1count = fig1count + 1;
%                                     legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
%                                 end
%                             end
                        elseif jointnum == 2
                            if abs(mean(dy(1:floor(size(dy,1)/2)))) > .85
                                if mean(dy(1:floor(size(dy,1)/2))) > 0
                                    figurecount = 1;
                                    fig1count = fig1count + 1;
                                    legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                else
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                end
                            else
                                figurecount = 3;
                                fig3count = fig3count + 1;
                                legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        elseif jointnum == 3
                            if transitions(1) < 500 && dd(1) < 0
                                figurecount = 3;
                                fig3count = fig3count + 1;
                                legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            elseif transitions(1) < 500 && dd(1) > 0
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                figurecount = 1;
                                fig1count = fig1count + 1;
                                legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        end
                    else
                        extflex = mean(dy);
                        [bb,~] = size(dy(dy>0));
                        [cc,~] = size(dy(dy<0));
                        imp = (dy(1)-mean(dy))/mean(dy);
                        if abs((bb-cc)/bb) < .6
                            figurecount = 3;
                            fig3count = fig3count + 1;
                            legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                        else
                            if extflex < 0 
                                figurecount = 1;
                                fig1count = fig1count + 1;
                                legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                figurecount = 2;
                                fig2count = fig2count + 1;
                                legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        end
                    end
                else
                    if fullstep
                        if jointnum == 1
                            if switches > 1 && abs(mean(dy(1:floor(size(dy,1)/2)))) < .9
                                figurecount = 3;
                                fig3count = fig3count + 1;
                                legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                if mean(dy(1:floor(size(dy,1)/2))) > 0
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                else
                                    figurecount = 1;
                                    fig1count = fig1count + 1;
                                    legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                end
                            end
                        elseif jointnum == 2
                            if abs(mean(dy(1:floor(size(dy,1)/2)))) > .85
                                if mean(dy(1:floor(size(dy,1)/2))) > 0
                                    figurecount = 1;
                                    fig1count = fig1count + 1;
                                    legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                else
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                                end
                            else
                                figurecount = 3;
                                fig3count = fig3count + 1;
                                legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        elseif jointnum == 3
                            if transitions(1) < 500 && dd(1) < 0
                                figurecount = 3;
                                fig3count = fig3count + 1;
                                legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            elseif transitions(1) < 500 && dd(1) > 0
                                    figurecount = 2;
                                    fig2count = fig2count + 1;
                                    legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                figurecount = 1;
                                fig1count = fig1count + 1;
                                legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        end
                    else
                        extflex = mean(dy);
                        [bb,~] = size(dy(dy>0));
                        [cc,~] = size(dy(dy<0));
                        imp = (dy(1)-mean(dy))/mean(dy);
                        if abs((bb-cc)/bb) < .6
                            figurecount = 3;
                            fig3count = fig3count + 1;
                            legender{1,figurecount}(fig3count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                        else
                            if extflex < 0 
                                figurecount = 1;
                                fig1count = fig1count + 1;
                                legender{1,figurecount}(fig1count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            else
                                figurecount = 2;
                                fig2count = fig2count + 1;
                                legender{1,figurecount}(fig2count,:) = cellstr(obj.musc_obj{relevant_muscles(i)}.muscle_name(4:end));
                            end
                        end
                    end
                end
                subplot(subplots,1,figurecount)
                if fullstep
                    plot(tnorm,mnorm)
                    if figurecount == 1
                        title(obj.joint_obj{jointnum}.name(4:end))
                    end
                    xlabel('Normalized step cycle')
                    ylabel('Muscle length (normalized)')
                else
                    plot(jointangle,mnorm)
                    xlabel(['Joint ',num2str(jointnum),' Angle (degrees)'])
                    ylabel('Muscle length (normalized)')
                end 
                legend(legender{1,figurecount},'Location','eastoutside')
                hold on
                extflex = 0;
            end
            if fullstep
                footer = [obj.joint_obj{jointnum}.name(4:end),'FullStep'];
            else
                footer = [obj.joint_obj{jointnum}.name(4:end),'JointRange'];
            end
            set(gcf,'Position',[500 50 900 900])
            axh = findobj( gcf, 'Type', 'Axes' );
            ax1 = get(axh(1),'Position');
            ax2= get(axh(2),'Position');
            ax3= get(axh(3),'Position');
            ax = [ax1;ax2;ax3];
            ax(:,3) = 1.5*min(ax(:,3));
            set(axh(1),'Position',ax(1,:))
            set(axh(2),'Position',ax(2,:))
            set(axh(3),'Position',ax(3,:))
            saveas(gcf,[pwd,'\OutputFigures\muscle_lengths_step\',footer,datestr(datetime('now'),'yyyymmdd'),'.png'])
        end
        %% Function: Plot Joint Passive Torque for a Step
        function plot_joint_passive_torque(obj)
            legendvec = cell(3,1);
            for i = 1:3
                subplot(2,1,1)
                [passive_joint_torque,passive_joint_motion] = compute_passive_muscle_torque(obj);
                legendvec{i,1} = obj.joint_obj{i}.name(4:end);
                xax = linspace(0,1,length(passive_joint_torque));
                plot(xax,passive_joint_torque(i,:))
                xlabel('Normalized Step Cycle')
                ylabel('Passive Joint Torque (Nm)')
                title('Passive Joint Torque - One Step')
                hold on
                subplot(2,1,2)
                plot(xax,passive_joint_motion(i,:))
                xlabel('Normalized Step Cycle')
                ylabel('Joint Motion (degrees)')
                hold on
            end
            subplot(2,1,1)
            legend(legendvec{:,1},'Location','eastoutside')
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            subplot(2,1,2)
            legend(legendvec{:,1},'Location','eastoutside')
            set(gcf,'Position',[200 200 1250 650])
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            saveas(gcf,[pwd,'\OutputFigures\Images\passive_joint_torque\','pjt',datestr(datetime('now'),'yyyymmdd'),'.png'])
        end
        %% Function: Optimize Forces for Torque Induction
        function [force,force2actuate,tau2,moment_output,fvals,telapsed,indiv_torques,exitflag] = optimize_forces(obj,toplot)
            tstart = tic;
            if nargin < 2
                toplot = 0;
            end
        
            xx = obj.sampling_vector;
            beg = obj.sampling_vector(1);
            ennd = obj.sampling_vector(end);
            
            tau2 = -obj.computer_active_joint_torque(0);
   
            nummuscles = size(obj.musc_obj,1);
                
            % Assemble an array Q, to act as a lower bound that ensures that Am values are possible
            Fmax = zeros(nummuscles,1);
            ub = zeros(nummuscles,length(tau2));
            lb = ub;
            fl = @(Lm,Lr,Lw) max(1-((Lm-Lr).^2./Lw^2),.7);
            m = length(obj.theta_motion(beg:ennd));
            n = length(obj.sampling_vector);
            dt = obj.dt_motion;
            for i=1:nummuscles
                muscle = obj.musc_obj{i};
                Fmax(i,1) = muscle.max_force;
                mL = muscle.muscle_length_profile;
                mV = muscle.muscle_velocity_profile;
                mVees(i,:) = mV;
                Lr = muscle.RestingLength;
                Lw = muscle.l_width;
                Al_musc = fl(mL,Lr,Lw);
                kp = muscle.Kpe;
                ks = muscle.Kse;
                b = muscle.damping;
                STmax = muscle.ST_max;
                if muscle.enabled
                    %lbHelper(i,:) = kp.*max(mL-Lr,0)+b.*mV;
                    lbHelper(i,:) = obj.passive_tension(:,i)';
                    ubHelper(i,:) = Al_musc.*STmax;
                    mTau = b/(ks+kp);
                    lbBnd(i,:) = [1/(1+dt/mTau),((ks*dt)/b)/(1+dt/mTau)]; % formulation found in notes 1-18-21, 4-15-21
                    ubBnd(i,:) = [1/(1+dt/mTau),((ks*dt)/b)/(1+dt/mTau),((ks*dt)/b)/(1+dt/mTau)]; % formulation found in notes 1-19-21, 4-15,21
                else
                    ub(i,:) = zeros(1,length(tau2));
                    lb(i,:) = ub(i,:);    
                end
            end
            x0 = zeros(size(ub,1),1);
            moment_output = zeros(nummuscles+2,size(xx,2),3);
            %%% REMEMBER TO REVERT SAMPLING VECTOR IN compute_muscle_moment_arm (line ~1355)
            for i=1:3
                moment_output(:,:,i) = compute_joint_moment_arms(obj,i,1);
            end
            force = zeros(nummuscles,length(xx));
            momentArmsHip = moment_output(:,:,1);
            momentArmsKnee = moment_output(:,:,2);
            momentArmsAnkle = moment_output(:,:,3);
            %options = optimoptions('fmincon','Display','none','Algorithm','sqp','OptimalityTolerance',1e-4,'Display','iter-detailed');
            options = optimoptions('fmincon','Display','none','Algorithm','sqp','OptimalityTolerance',1e-5);
            exitflag = -10.*ones(1,length(xx));
            fvals = -1.*ones(3,length(xx));
            
            ub = Inf(38,1); lbMat = zeros(38,length(xx)); count = 1;
            lb = zeros(38,1);
            for j = 1:length(xx)
                Aeq = [momentArmsHip(1:nummuscles,j),momentArmsKnee(1:nummuscles,j),momentArmsAnkle(1:nummuscles,j)]';
                beq = 1000.*tau2(j,:);
                if j>1
                    x0 = force(:,j-1);
                    fun = @(x) sum((x./Fmax).^2)+0*sum((x-x0).^2);
                else
                    fun = @(x) sum((x./Fmax).^2);
                end
                lb = max(min(lbBnd(:,1).*x0+lbBnd(:,2).*lbHelper(:,j),Fmax),obj.passive_tension(j,:)');
                %ub = max(min(ubBnd(:,1).*x0+ubBnd(:,3).*ubHelper(:,j)+ubBnd(:,3).*lbHelper(:,j),Fmax),0);
                lbLog(j,:) = lbBnd(:,1).*x0+lbBnd(:,2).*lbHelper(:,j);
                ubLog(j,:) = ubBnd(:,1).*x0+ubBnd(:,3).*ubHelper(:,j)+ubBnd(:,3).*lbHelper(:,j);
                if any(lb > ub)
                    lb(lb>ub) = .9999*ub(lb>ub);
                end
                [force(: ,j),~,exitflag(j)] = fmincon(fun,x0,[],[],Aeq,beq,lb,ub,[],options);
                fvals(:,j) = Aeq*force(:,j)-beq';
                temp(:,j) = [sum((force(:,j)./Fmax).^2),sum((force(:,j)-x0).^2)];
                lbMat(:,j) = lb;
                %bb = 1;
                %disp(num2str(count))
                count = count +1;
            end
            
            muscZones = zoning_sorter(obj.original_text,6);
            bb = max((force./Fmax),[],2);
            muscZones(:,3) = num2cell(bb);
            cc = sortrows(muscZones,2);
            
            indiv_torques = zeros(nummuscles,length(xx),3);
            for ii = 1:3
                indiv_torques(1:nummuscles,:,ii) = moment_output(1:nummuscles,:,ii)./1000.*force;
            end
            
            force = force';
            force2actuate = force-obj.passive_tension;
%             % For outputting time
            telapsed = toc(tstart);
            if telapsed>60
                min_time = num2str(floor(telapsed/60));
                sec = num2str(round(mod(telapsed,60)));
            else
                min_time = num2str(0);
                sec = num2str(round(telapsed,2));
            end
%             disp(['Force Optimization Time:',' (',min,'m ',sec,'s)'])
            %%
            if toplot
                forces2plot = [];
                legender = {};
                for ii = 1:size(force,1)
                    if sum(force(ii,:)) == 0
                    else
                        legender{1,end+1} = obj.musc_obj{ii}.muscle_name(4:end);
                        forces2plot = [forces2plot;force(ii,:)];
                    end
                end

                forcefig = figure;
                    set(forcefig,'Position',[200 100 1500 900])
                    forcesub = subplot(4,1,1);
                        CM = hsv(size(forces2plot,1));
                        for jj = 1:size(forces2plot,1)
                            plot(forces2plot(jj,:)','Color',CM(jj,:),'LineWidth',2)
                            hold on
                        end
                        %title(mintypes{optset})
                        legend(legender,'Location','eastoutside')
                    colorforces2 = subplot(4,1,2);
                        CM = lines(size(forces2plot,1));
                        for jj = 1:size(forces2plot,1)
                            plot(forces2plot(jj,:)','Color',CM(jj,:),'LineWidth',2)
                            hold on
                        end
                        %colorforces2.Position = [colorforces2.Position(1) colorforces2.Position(2) forcesub.Position(3) forcesub.Position(4)];
                        legend(legender,'Location','westoutside')
                        drawnow
                        colorforces2.Position = [forcesub.Position(1) colorforces2.Position(2) forcesub.Position(3) forcesub.Position(4)];
                    jointsub = subplot(4,1,3);
                        plot(obj.theta_motion(xx,:),'LineWidth',2)
                        legend({'Hip','Knee','Ankle'},'Location','eastoutside')
                        drawnow;
                        jointsub.Position = [jointsub.Position(1) jointsub.Position(2) forcesub.Position(3) forcesub.Position(4)];
                    torquesub = subplot(4,1,4);
                        plot(tau2,'LineWidth',2)
                        legend({'Hip','Knee','Ankle'},'Location','eastoutside')
                        drawnow;
                        torquesub.Position = [torquesub.Position(1) torquesub.Position(2) forcesub.Position(3) forcesub.Position(4)];
                    clear CM

                easyforce = figure;
                    set(easyforce,'Position',[200 100 1500 900])
                    numforces = size(legender,2);
                    for ii = 1:4
                        subplot(4,1,ii)
                        if 5*ii > length(legender)
                            forceindex = ((1+5*(ii-1)):numforces);
                        else
                            forceindex = (1+5*(ii-1)):(5*ii);
                        end
                        plot(forces2plot(forceindex,:)','LineWidth',2.5)
                        hold on
                        ylimit = ylim;
                        ylimit = ylimit(2);
                        shadeheight = ones(1,100)*ylimit;
                        shadelength  = 1:100;
                        area(shadelength(1:50),shadeheight(1:50),'FaceAlpha',.05,'EdgeAlpha',0,'FaceColor',[1 0 0])
                        area(shadelength(50:100),shadeheight(50:100),'FaceAlpha',.05,'EdgeAlpha',0,'FaceColor',[0 0 1])
                        legend(legender(forceindex),'Location','eastoutside','FontSize',12);
                        xlabel('Percent Stride (%)')
                        ylabel('Force (N)')
                        grid on
                        if ii == 1
                            firstsubposition = get(gca,'Position');
                            title('Optimized')
                        else
                            currentsubposition = get(gca,'Position');
                            currentsubposition = [currentsubposition(1:2),firstsubposition(3:4)];
                            set(gca,'Position',currentsubposition);
                        end
                    end
                    clear numforces ylimit shadeheight shadelength firstsubposition currentsubposition
            end
        end    
        %% Function: Forward Dynamics: Compute mass matrix for given limb configuration
        function [outM,Mprime] = compute_mass_matrix(obj,theta)
                % Calculates the inertia matrix for a mulit-link actuator
                % From Murray, Li, Sastry 94, p.176 (4.29)
                
                % Adjoint transformation, p.55, (2.58)
                Adj_func = @(eX) [[eX(1:3,1:3) wedge(eX(1:3,4))*eX(1:3,1:3)];[zeros(3,3) eX(1:3,1:3)]];
                % Exponential mapping function. Definition on p.42 (2.36)
                exp_func = @(R,w,v,th) [[R (eye(3)-R)*(cross(w,v))+w*w'*v*th];[zeros(1,3) 1]];
                
                gWH = [[obj.body_obj{1}.CR [0;0;0]];[zeros(1,3) 1]];
                gHF = [[obj.body_obj{2}.CR obj.body_obj{2}.position'];[zeros(1,3) 1]];
                gFT = [[obj.body_obj{3}.CR obj.body_obj{3}.position'];[zeros(1,3) 1]];
                gTX = [[obj.body_obj{4}.CR obj.body_obj{4}.position'];[zeros(1,3) 1]];

                femCOMw = gWH*gHF*[obj.body_obj{2}.com';1];
                tibCOMw = gWH*gHF*gFT*[obj.body_obj{3}.com';1];
                fotCOMw = gWH*gHF*gFT*gTX*[obj.body_obj{4}.com';1];

                gslfem = [[obj.body_obj{2}.Cabs femCOMw(1:3)];[zeros(1,3) 1]];
                gsltib = [[obj.body_obj{3}.Cabs tibCOMw(1:3)];[zeros(1,3) 1]];
                gslfot = [[obj.body_obj{4}.Cabs fotCOMw(1:3)];[zeros(1,3) 1]];
                
                m1 = obj.body_obj{2}.mass/1000;
                m2 = obj.body_obj{3}.mass/1000;
                m3 = obj.body_obj{4}.mass/1000;

                % Computing the moments of inertia
                % These are just approximations, assuming femur and tibia are solid rods and the foot is a box
                % Measurements for dimensions approximated by inspection using Animatlab
                MOI = zeros(3,3,3);
                % Femur 
                L1 = obj.body_obj{2}.length; r1 = 0.049654*L1; MOI(1,1,1) = (1/12)*m1*L1^2; MOI(2,2,1) = (1/12)*m1*L1^2; MOI(3,3,1) = (1/2)*m1*r1^2;
                % Tibia
                L2 = obj.body_obj{3}.length; r2 = 0.03061*L2; MOI(1,1,2) = (1/12)*m2*L2^2; MOI(2,2,2) = (1/12)*m2*L2^2; MOI(3,3,2) = (1/12)*m2*L2^2;
                % Tibia
                d = obj.body_obj{4}.length; h = 0.13473*d; w = 0.44908*d; MOI(1,1,3) = (1/12)*m3*(d^2+h^2); MOI(2,2,3) = (1/12)*m3*(d^2+w^2); MOI(3,3,3) = (1/12)*m3*(h^2+w^2);
                
                % Calculate the general inertia matrices, p.172 (bottom)
                gifem = [[eye(3).*[m1;m1;m1] zeros(3,3)];[zeros(3,3) eye(3).*[MOI(1,1,1);MOI(2,2,1);MOI(3,3,1)]]];
                gitib = [[eye(3).*[m2;m2;m2] zeros(3,3)];[zeros(3,3) eye(3).*[MOI(1,1,2);MOI(2,2,2);MOI(3,3,2)]]];
                gifot = [[eye(3).*[m3;m3;m3] zeros(3,3)];[zeros(3,3) eye(3).*[MOI(1,1,3);MOI(2,2,3);MOI(3,3,3)]]];
                
                Mprime(:,:,1) = inv(Adj_func(gslfem))'*gifem*inv(Adj_func(gslfem));
                Mprime(:,:,2) = inv(Adj_func(gsltib))'*gitib*inv(Adj_func(gsltib));
                Mprime(:,:,3) = inv(Adj_func(gslfot))'*gifot*inv(Adj_func(gslfot));
                
                axesMat = joint_axes_from_angles(obj,theta);
                jointMat = joint_pos_on_demand(obj,theta);
                
                Jac = [-cross(axesMat,jointMat);axesMat];

                e1 = exp_func(obj.axis_angle_rotation(theta(1),Jac(4:6,1)),Jac(4:6,1),Jac(1:3,1),theta(1));
                e2 = exp_func(obj.axis_angle_rotation(theta(2),Jac(4:6,2)),Jac(4:6,2),Jac(1:3,2),theta(2));
                e3 = exp_func(obj.axis_angle_rotation(theta(3),Jac(4:6,3)),Jac(4:6,3),Jac(1:3,3),theta(3));
                
                J1 = [inv(Adj_func(e1*gslfem))*Jac(:,1) zeros(6,1) zeros(6,1)];
                J2 = [inv(Adj_func(e1*e2*gsltib))*Jac(:,1) inv(Adj_func(e2*gsltib))*Jac(:,2) zeros(6,1)];
                J3 = [inv(Adj_func(e1*e2*e3*gslfot))*Jac(:,1) inv(Adj_func(e2*e3*gslfot))*Jac(:,2) inv(Adj_func(e3*gslfot))*Jac(:,3)];

                outM = J1'*gifem*J1+J2'*gitib*J2+J3'*gifot*J3;
        end
        %% Function: Forward Dynamics: Compute gravitational force vector for a given limb configuration
        function outG = compute_gravity_vector(obj,theta)
            
            scale = 1000;

            jointv = joint_axes_from_angles(obj,theta)*10;
            jointMat = joint_pos_on_demand(obj,theta);
            COMmat = [com_pos_on_demand(obj,theta,2) com_pos_on_demand(obj,theta,3) com_pos_on_demand(obj,theta,4)];
            
            % This is an array of points defining which bodies to take the moment arm for
            % For example, the first page is the femur COM about the hip and then the knee about the hip
            armPoints(:,:,1) = [COMmat(:,1) jointMat(:,1) jointMat(:,2) jointMat(:,1)];
            armPoints(:,:,2) = [COMmat(:,2) jointMat(:,2) jointMat(:,3) jointMat(:,2)];
            armPoints(:,:,3) = [COMmat(:,3) jointMat(:,3) zeros(3,1) zeros(3,1)];
            
            for jointNum = 1:3
                count = 1;
                for armNum = 1:2:3
                    % For each pair of points in the armPoints array, find the moment arm
                    pointRel = zeros(2,3);
                    pointProjection = pointRel;
                    pointRel = [armPoints(:,armNum,jointNum) armPoints(:,armNum,jointNum)-[0 ;-.05; 0]]'-armPoints(:,armNum+1,jointNum)';
                    pointProjection(1,:) = (pointRel(1,:) - (dot(pointRel(1,:),jointv(:,jointNum)')/norm(jointv(:,jointNum)')^2)*jointv(:,jointNum)')'+jointMat(:,jointNum);
                    pointProjection(2,:) = (pointRel(2,:) - (dot(pointRel(2,:),jointv(:,jointNum)')/norm(jointv(:,jointNum)')^2)*jointv(:,jointNum)')'+jointMat(:,jointNum);
                    pointProjection = scale*pointProjection;
                    %Matrix of muscle projections perpendicular to the joint axis
                        fflatmat = pointProjection(2,:)-pointProjection(1,:);
                    %Storing the moment arms in a matrix
                        momentArmLong = cross(fflatmat,jointv(:,jointNum));
                        scaledjoint = scale*jointMat(:,jointNum)';
                        PA2 = [pointProjection(1,:);scaledjoint];
                        PB2 = [pointProjection(2,:);scaledjoint+momentArmLong];
                    %lineINtersect3D gives the scaled vector of the moment arms for each segment
                        momentArm = lineIntersect3D(obj,PA2,PB2);
                    %A modifier to determine whether the moment arm is positive or negative
                        sig_momentarm = momentArm-scaledjoint;
                        sig_muscleprojection = (pointProjection(1,:)-pointProjection(2,:));
                        signal2 = sign(dot(cross(sig_momentarm,sig_muscleprojection),jointv(:,jointNum)'));
                        momentMat(jointNum,count) = signal2*norm(sig_momentarm)/1000;
                        count = count +1;
                end
            end

            m = zeros(1,4);
            for ii = 1:length(obj.body_obj)
                % Body mass in kg
                m(ii) = obj.body_obj{ii}.mass./1000;
            end
                
            outG(1) = m(2)*momentMat(1,1)+(m(3)+m(4))*momentMat(1,2);
            outG(2) = m(3)*momentMat(2,1)+m(4)*momentMat(2,2);
            outG(3) = m(4)*momentMat(3,1);
            
            outG = outG*9.8;
        end
        %% Function: Forward Dynamics: Compute coriolis matrix
        function [outC,outM] = compute_coriolis_matrix(obj,theta,thetaDot)
            % Make the Coriolis term according to MLS 1994 p. 176 (4.30)
            % Input: A cell(3,3): The adjoint Transformation according to MLS 1994 p. 176 (4.27)
            % Input: Jac double(6,3): spatial manipulator Jacobian. Each column is the joint twist
            % Input: M double(6,6): generalized inertia matrix following form MLS 1994 p.176 (4.28)
            % Input: theta_dot double(3,1): velocity of joint angle input
        
            A = compute_adjoint_transformation_matrix(obj,theta);
            axesMat = joint_axes_from_angles(obj,theta);
            jointMat = joint_pos_on_demand(obj,theta);

            Jac = [-cross(axesMat,jointMat);axesMat];
            
            [outM,M] = compute_mass_matrix(obj,theta);
            
            outC = zeros(3,3);
            for i = 1:3
                for j = 1:3
                    temp = 0;
                    for k = 1:3
                        dMdT_1 = 0;
                        for l = max(i,j):3
                            dMdT_1 = dMdT_1+bracket(A{k,i}*Jac(:,i),Jac(:,k))'*A{l,k}'*M(:,:,l)*A{l,j}*Jac(:,j)+Jac(:,i)'*A{l,i}'*M(:,:,l)*A{l,k}*bracket(A{k,j}*Jac(:,j),Jac(:,k));
                        end
                        dMdT_2 = 0;
                        for l = max(i,k):3
                            dMdT_2 = dMdT_2+bracket(A{j,i}*Jac(:,i),Jac(:,j))'*A{l,j}'*M(:,:,l)*A{l,k}*Jac(:,k)+Jac(:,i)'*A{l,i}'*M(:,:,l)*A{l,j}*bracket(A{j,k}*Jac(:,k),Jac(:,j));
                        end
                        dMdT_3 = 0;
                        for l = max(k,j):3
                            dMdT_3 = dMdT_3+bracket(A{i,k}*Jac(:,k),Jac(:,i))'*A{l,i}'*M(:,:,l)*A{l,j}*Jac(:,j)+Jac(:,k)'*A{l,k}'*M(:,:,l)*A{l,i}*bracket(A{i,j}*Jac(:,j),Jac(:,i));
                        end
                        temp = temp+(1/2)*(dMdT_1+dMdT_2-dMdT_3)*thetaDot(k);
                    end
                    outC(i,j) = temp;
                end
            end
        end
        %% Function: Forward Dynamics: Compute Adjoint Transformation Matrix
        function outA = compute_adjoint_transformation_matrix(obj,theta)
            % Create the adjoint transformation matrix according to MLS 1994 p. 176 (4.27)
            Adj_func = @(eX) [[eX(1:3,1:3) wedge(eX(1:3,4))*eX(1:3,1:3)];[zeros(3,3) eX(1:3,1:3)]];
            exp_func = @(R,w,v,th) [[R (eye(3)-R)*(cross(w,v))+w*w'*v*th];[zeros(1,3) 1]];
            
            axesMat = joint_axes_from_angles(obj,theta);
            jointMat = joint_pos_on_demand(obj,theta);

            Jac = [-cross(axesMat,jointMat);axesMat];
            
            for hh = 1:3
                for kk = 1:3
                    if hh>kk
                        expMat = eye(4);
                        for tt = kk+1:hh
                            temp1  = exp_func(obj.axis_angle_rotation(theta(tt),Jac(4:6,tt)),Jac(4:6,tt),Jac(1:3,tt),theta(tt));
                            expMat = expMat*temp1; 
                        end
                        outA{hh,kk} = inv(Adj_func(expMat));
                    elseif hh==kk
                        outA{hh,kk} = eye(6);
                    elseif hh<kk
                        outA{hh,kk} = zeros(6,6);
                    end
                end
            end
        end
        %% Function: Pedotti Optimization
        function results_cell = pedotti_optimization(obj)
            results_cell = cell(10,2);
            results_cell(:,1) = {'';
                                 'Total Forces Needed';
                                 'Force That the MN Needs to Generate';
                                 'Torques';
                                 'Moment Arm output';
                                 'Fval';
                                 'Opt Time';
                                 'Indiv Torques';
                                 'Exitflags';
                                 'Muscles Ranked by Force'};
            results_cell(1,2) = {' '};
            num_muscles = length(obj.musc_obj);
            for i=1:num_muscles
                muscle = obj.musc_obj{i};
                Fmax(i,1) = muscle.max_force;
            end

            fun = @(x) sum((x./Fmax).^2);

            [forces,force2actuate,tau,moment_output,fvals,telapsed,indiv_torques,exitflags] = obj.optimize_forces(0);
            
            for ii = 1:38
                areaunder{ii,1} = obj.musc_obj{ii}.muscle_name(4:end);
                areaunder{ii,2} = trapz(forces(:,ii));
            end
            
            results_cell{2,2} = forces;
            results_cell{3,2} = force2actuate;
            results_cell{4,2} = tau;
            results_cell{5,2} = moment_output;
            results_cell{6,2} = fvals';
            results_cell{7,2} = telapsed;
            results_cell{8,2} = indiv_torques;
            results_cell{9,2} = exitflags';
            results_cell{10,2} = sortrows(areaunder,2,'descend');
        end
        %% Function: Visualize Force Activation over Stride
        function vis_walking(obj,forces,to_save)
            % All forces must be positive for this to work. Line widths are based on activation level, can't be negative.
            %zone_colors = [0,255,0,255,0,0,235,210,52,0,0,255,255,0,255,0,255,255;217,255,217,255,215,215,235,231,203,234,234,255,255,222,255,212,255,255]./255;
            zone_colors = [0,255,0,255,0,0,235,210,52,0,0,255,255,0,255,0,255,255;255.*ones(1,18)]./255;
            zc_inds = 1:3:length(zone_colors);   
            for ii = 1:6
                zcs(:,:,ii) = interp1([100;0],zone_colors(:,zc_inds(ii):zc_inds(ii)+2),linspace(100,0,155),'pchip');
            end
            muscZones = zoning_sorter(obj.original_text,6);
            ind = 117:3039;
            % Create the color map. Need to trim it down a bit to be visible on the ends of the color range
            tempcolor = hot(256);
            trim = floor(.2*length(tempcolor));
            tempcolor = tempcolor(trim:(end-trim),:);
            CM = flipud(tempcolor);

            for musnum = 1:size(forces,2)
                force_vec = (forces(ind,musnum) - min(forces(ind,musnum)))/(max(forces(ind,musnum)) - min(forces(ind,musnum)));
                force_vec(isnan(force_vec)) = 0;
                widths(:,musnum) = (6-.1).*force_vec+.1;
                colors(:,musnum) = 1+max(floor((length(CM)-1)*force_vec),0);
            end

            aa = figure('Position',[800,266,998,684]);
            colormap(aa,CM);
            num_muscles = length(obj.musc_obj);
            savePath = 'G:\My Drive\Rat\MeetingFiles\Meeting_20210503\LegAnim';
            sPos = [[0.0666,0.108,0.524,0.816];[0.563,0.11,0.335,0.815]];
            filePath = tempname(savePath);
            for timecount = ind(1):50:ind(end)    
                for splot = 1:2
                    if timecount == ind(1)
                        subplots{splot} = subplot(1,2,splot);
                    else
                        cla(subplots{splot})
                    end
                    for musnum = 1:num_muscles
                        muscle = obj.musc_obj{musnum};
                        musclemat =zeros(size(muscle.pos_attachments,1),3);
                        for hh = 1:size(muscle.pos_attachments,1)
                            temp = 1000*muscle.pos_attachments{hh,4}(timecount,:);
                            musclemat(hh,:) = temp;
                        end
                         plot3(subplots{splot},musclemat(:,1),musclemat(:,2),musclemat(:,3),'Color',CM(colors(timecount-ind(1)+1,musnum),:),'Linewidth',widths(timecount-ind(1)+1,musnum))
%                            plot3(subplots{splot},musclemat(:,1),musclemat(:,2),musclemat(:,3),'Color',zcs(colors(timecount-ind(1)+1,musnum),:,muscZones{musnum,2}),'Linewidth',widths(timecount-ind(1)+1,musnum))
                        hold on
                    end
                    jointmat = zeros(4,3);
                    for jointnum = 1:3
                        jointmat(jointnum,:) = 1000*obj.joint_obj{jointnum}.sim_position_profile(timecount,:);
                    end
                    jointmat(4,:) = 1000*obj.musc_obj{20}.pos_attachments{end,4}(timecount,:);
                    plot3(subplots{splot},jointmat(:,1),jointmat(:,2),jointmat(:,3),'m--','LineWidth',2,'Marker','o','MarkerFaceColor','k','MarkerSize',5)
                     grid on
                        % 1x
                        xlim([-40 35]);
                        ylim([-65 15]);
                        zlim([0 30]);
                        %10x
%                         xlim([-400 304.08]);
%                         ylim([-604.08 150]);
%                         zlim([0 300]);
                        % 50x
%                         xlim([-2000 1500]);
%                         ylim([-3000 750]);
%                         zlim([0 1500]);
                    pbaspect([sum(abs(get(gca,'XLim'))) sum(abs(get(gca,'YLim'))) sum(abs(get(gca,'ZLim')))]/norm([sum(abs(get(gca,'XLim'))) sum(abs(get(gca,'YLim'))) sum(abs(get(gca,'ZLim')))]))
                    if timecount == ind(1)
                        switch splot
                            case 1
                                view([0 90]); 
                            case 2
                                view([90 0]);
                                camroll(90)
                                cbar = colorbar('Ticks',[0;1],'TickLabels',{'low act','high act'},'FontSize',15,'Position',[0.881,0.109,0.0199,0.815]);
                        end
                    end
                    drawnow
                end
                if timecount == ind(1)
                    set(subplots{1},'Position',[0.0666,0.108,0.524,0.816])
                    set(subplots{2},'Position',[0.563,0.11,0.335,0.815])
                end
%                 set(cbar,'Position',[0.881,0.109,0.0199,0.815])
%                 pause(.001)
               % saveas(aa,[savePath,zoneNames{zoneNum},'_',num2str(lengthNum),'.png'])
                % Capture the plot as an image
                if to_save
                    frame = getframe(aa); 
                    im = frame2im(frame); 
                    [imind,cm] = rgb2ind(im,256); 
                    % Write to the GIF File 
                    if timecount == ind(1) 
                      imwrite(imind,cm,[filePath,'.gif'],'gif', 'Loopcount',inf,'DelayTime',.54e-3); 
                    else 
                      imwrite(imind,cm,[filePath,'.gif'],'gif','WriteMode','append','DelayTime',.54e-3); 
                    end 
                end
                hold off
            end
        end
    end
end