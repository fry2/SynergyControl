function scaledPathName = legScaler(inPath,factor,inVals)
    % For an input ASIM or APROJ path, scale the leg by a certain factor. Writes a new file, located at scaledPathName, based on the contents of the file at inPath
    % Input: inPath (char): full path to the project to be scaled
    % Input: factor (double): value for scaling (e.g. 2 = 2x larger)
    
    % Import project text
    if ~ischar(inPath)
        inPath = char(inPath);
    end
    
    if contains(inPath,'.aproj')
        is_proj = 1;
    elseif contains(inPath,'.asim')
        is_proj = 0;
    else
        error('inPath doesn''t point to an .asim or .aproj file')
    end
        
    if nargin < 3
        inVals = [];
    end
    
    inText = importdata(inPath);
    inText = changeWalkingWaveforms(inText,factor,is_proj);
    % Find organism position
    if is_proj
        tempind = find(contains(inText,'<Organism>'))+17;
        counter = 1; orgpos = zeros(1,3);
        for ii = tempind+1:tempind+3
            orgpos(counter) = double(extractBetween(string(inText{ii}),'Actual="','"'));
            counter = counter + 1;
        end
    else
        tempind = find(contains(inText,'<Organism>'))+15;
        orgpos = cellfun(@str2num,extractBetween(inText{tempind},'"','"'))';
    end
    
    % Update bone masses
    boneNames = {'Pelvis','LH_Femur','LH_Tibia','LH_Foot'};
    massScale = char(extractBetween(string(inText{contains(inText,'MassUnits')}),'>','</'));
    switch massScale
        case 'Grams'
            ms = 1e-3;
        case 'Kilograms'
            ms = 1;
        otherwise
            error('Check the mass scale for the project, only prepared to accept Grams or Kilograms');
    end
    for ii = 1:length(boneNames)
        if is_proj
            boneInd = find(contains(inText,['<Name>',boneNames{ii},'</Name>']));
            massInd = find(contains(inText(boneInd:end),'<Mass'),1,'first')+boneInd-1;
            mass = double(extractBetween(string(inText{massInd}),'Actual="','"')); 
            inText{massInd} = replacePositionInLine(inText{massInd},factor^3*mass);
        else
%             disp('CHECK THIS, YOU CHANGED THE MASSES AFTER DOING THIS.')
%             disp('PEL 31.824G, FEM 14.141G, TIB 3.342G, FOOT 1.571G')
%             keyboard
            boneInd = find(contains(inText,['<Name>',boneNames{ii},'</Name>']));
            massInd = find(contains(inText(boneInd:end),'<Mass>'),1,'first')+boneInd-1;
            % In order to make the model work at different scales, we need to scale the actual bone meshes
            % Animatlab refuses to acknowledge this when it's done programmtically. As a workaround, we can shrink those meshes to a tiny scale
            % so they won't collide with each other and then perform all the usual kinematics with the muscles and joints. The one catch is that
            % Animatlab automatically calculates new bone masses based on this super-shrunk size. So we need to scale the bone masses back to their original masses
            % before we actually do this "new" scaling in the rest of the program
            boneNameTrim = regexprep(boneNames{ii},'_','','emptymatch');
            boneFactor = double(extractBetween(string(inText{find(contains(inText,[boneNameTrim,'_Convex.osg</ConvexMeshFile>']))+1}),'x="','"'));
            mass = double(extractBetween(string(inText{massInd}),'>','</'));
            %actualMass = mass/(boneFactor^3);
            actualMass = mass;
            inText{massInd} = ['<Mass>',num2str(factor^3*actualMass),'</Mass>'];
        end
    end
    
    % Update all scaling factors in the project text
    if is_proj
        scaleIndVec = find(contains(inText,'<Scale>'));
        for ii = 1:length(scaleIndVec)
            scaleInd = scaleIndVec(ii);
            for jj = scaleInd+1:scaleInd+3
                 inText{jj} = replacePositionInLine(inText{jj},factor);
            end
        end
    else
        scaleIndVec = find(contains(inText,'<Scale'));
        for ii = 1:length(scaleIndVec)
            inText{scaleIndVec(ii)} = replaceBetween(inText{scaleIndVec(ii)},'"','"',num2str(factor));
        end
    end
    
    % Update all local positions based on new scale
    if is_proj
        rbIndVec = find(contains(inText,'<RigidBody>'));
        for ii = 1:length(rbIndVec)
            rb1 = rbIndVec(ii);
            typeInd = find(contains(inText(rb1:end),'<Type>'),1,'first')+rb1-1;
            bodyID = contains({'<Name>Pelvis</Name>';'<Name>LH_Femur</Name>';'<Name>LH_Tibia</Name>';'<Name>LH_Foot</Name>'},inText{rb1+1});
            if ~contains(inText{rb1+3},'<Type>LinearHillMuscle</Type>') % If you're not a muscle (attachment, bone, etc)
                % Update Local position
                lpInd = find(contains(inText(rb1:end),'<LocalPosition>'),1,'first')+rb1-1;
                counter = 1; newlp = zeros(1,3);
                for jj = lpInd+1:lpInd+3
                    oldVal = double(extractBetween(string(inText{jj}),'Actual="','"'));
                    inText{jj} = replacePositionInLine(inText{jj},factor*oldVal);
                    newlp(counter) = factor*oldVal;
                    counter = counter + 1;
                end
                % Update COMs
                if contains(inText{typeInd},'Mesh')
                    comInd = find(contains(inText(rb1:end),'<COM>'),1,'first')+rb1-1;
                    for jj = comInd+1:comInd+3
                        oldVal = double(extractBetween(string(inText{jj}),'Actual="','"'));
                        inText{jj} = replacePositionInLine(inText{jj},factor*oldVal);
                    end
                end
                % Update Local Matrix
                lpMatInd = find(contains(inText(rb1:end),'<LocalMatrix'),1,'first')+rb1-1;
                oldMat = sscanf(char(extractBetween(string(inText{lpMatInd}),'>','</')), '%g,')';
                oldMat(13:15) = 10*newlp;
                newMat = [];
                for jj = 1:length(oldMat)
                    if oldMat(jj) == 0
                        newMat = [newMat,'0,'];
                    elseif oldMat(jj) == 1
                        newMat = [newMat,'1,'];
                    elseif oldMat(jj) == -1
                        newMat = [newMat,'-1,'];
                    else
                        newMat = [newMat,[sprintf('%.6f',oldMat(jj)),',']];
                    end
                end
                %newMat = sprintf('%.6f,',oldMat);
                newMat = newMat(1:end-1);
                inText{lpMatInd} = replaceBetween(inText{lpMatInd},'>','</',newMat);
                % Update Dragger radius
                dragInd = find(contains(inText(rb1:end),'<DraggerSize'),1,'first')+rb1-1;
                dragSize = double(extractBetween(string(inText{dragInd}),'Actual="','"'));
                if dragSize ~= -1
                    inText{dragInd} = replacePositionInLine(inText{dragInd},factor*dragSize);
                end
                % Update Attachment Radius
                if contains(inText{typeInd},'Attachment')
                    radiusInd = find(contains(inText(rb1:end),'<Radius'),1,'first')+rb1-1;
                    radiusSize = double(extractBetween(string(inText{radiusInd}),'Actual="','"'));
                    inText{radiusInd} = replacePositionInLine(inText{radiusInd},factor*radiusSize);
                end
            else % is_proj, it is a muscle
                % Change the muscle parameters
                    ksInd = find(contains(inText(rb1:end),'Kse'),1,'first')+rb1-1;
                    kpInd = find(contains(inText(rb1:end),'Kpe'),1,'first')+rb1-1;
                    bInd = kpInd+1; % This is easier than searching for the second B Value (since the ST curve uses a B Value, too)
                    muscParamInds = [ksInd,kpInd,bInd];
                    for jj = 1:length(muscParamInds)
                        oldVal = double(extractBetween(string(inText{muscParamInds(jj)}),'Actual="','"'));
                        inText{muscParamInds(jj)} = replacePositionInLine(inText{muscParamInds(jj)},factor*oldVal);
                        newVEs(jj) = factor*oldVal;
                    end
                % For individual Fmax scaling
                    fmInd = find(contains(inText(rb1:end),'MaximumTension'),1,'first')+rb1-1;
                    oldFM = double(extractBetween(string(inText{fmInd}),'Actual="','"'));
                    if factor < 1
                        newFM = factor^2*oldFM;
                    else
                        newFM = factor^3*oldFM;
                    end
                    inText{fmInd} = replacePositionInLine(inText{fmInd},newFM);
                % Change LT params
                    ltInd = find(contains(inText(rb1:end),'<LengthTension>'),1,'first')+rb1-1;
                    lrInd = find(contains(inText(ltInd:end),'RestingLength'),1,'first')+ltInd-1; % Resting Length
                    lwInd = find(contains(inText(ltInd:end),'Lwidth'),1,'first')+ltInd-1; % Lwidth
                    llInd = find(contains(inText(ltInd:end),'LowerLimitScale'),1,'first')+ltInd-1; % Lower Limit
                    ulInd = find(contains(inText(ltInd:end),'UpperLimitScale'),1,'first')+ltInd-1; % Upper Limit
                    ltIndVec = [lrInd, lwInd, llInd, ulInd];
                    for jj = 1:length(ltIndVec)
                        oldVal = double(extractBetween(string(inText{ltIndVec(jj)}),'Actual="','"'));
                        inText{ltIndVec(jj)} = replacePositionInLine(inText{ltIndVec(jj)},factor*oldVal);
                    end
                % Change ST params
                    newSTmax = (1+newVEs(2)/newVEs(1))*newFM;
                    newYoff = -.01*newSTmax;
                    stInd = find(contains(inText(rb1:end),'<StimulusTension>'),1,'first')+rb1-1;
                    stMaxInd = find(contains(inText(stInd:end),'<B Value'),1,'first')+stInd-1; % STmax
                    yOffInd = find(contains(inText(stInd:end),'<D Value'),1,'first')+stInd-1; % Lwidth
                    inText{stMaxInd} = replacePositionInLine(inText{stMaxInd},newSTmax);
                    inText{stMaxInd-4} = replacePositionInLine(inText{stMaxInd-4},newSTmax);
                    inText{yOffInd} = replacePositionInLine(inText{yOffInd},newYoff);
            end
        end
        jIndVec = find(contains(inText,'<Joint>'));
        % Update JointSize
        for ii = 1:length(jIndVec)
            j1 = jIndVec(ii);
            jSizeInd = find(contains(inText(j1:end),'<Size'),1,'first')+j1-1;
            jSize = double(extractBetween(string(inText{jSizeInd}),'Actual="','"'));
            inText{jSizeInd} = replacePositionInLine(inText{jSizeInd},factor*jSize);
        end
    else % Is simulation file
        rbIndVec = find(contains(inText,'<RigidBody>'));
        for ii = 1:length(rbIndVec) % For each rigid body
            rb1 = rbIndVec(ii);
            typeInd = find(contains(inText(rb1:end),'<Type>'),1,'first')+rb1-1;
            if ~contains(inText{rb1+3},'<Type>LinearHillMuscle</Type>') % is_sim, obj is not a muscle (bone, attachment, etc)
                % Update Local Positions
                lpInd = find(contains(inText(rb1:end),'<Position'),1,'first')+rb1-1;
                oldVals = cellfun(@str2num,extractBetween(inText{lpInd},'"','"'))';
                newVals = factor.*oldVals;
                for jj = 1:3
                    quoteLocs = strfind(inText{lpInd},'"');
                    qts = quoteLocs(2*jj-1:2*jj);
                    inText{lpInd} = replaceBetween(inText{lpInd},qts(1)+1,qts(2)-1,num2str(newVals(jj)));
                end
                % Update COMs
                if contains(inText{typeInd},'Mesh')
                    comInd = find(contains(inText(rb1:end),'<COM'),1,'first')+rb1-1;
                    oldVals = cellfun(@str2num,extractBetween(inText{comInd},'"','"'))';
                    newVals = factor.*oldVals;
                    for jj = 1:3
                        quoteLocs = strfind(inText{comInd},'"');
                        qts = quoteLocs(2*jj-1:2*jj);
                        inText{comInd} = replaceBetween(inText{comInd},qts(1)+1,qts(2)-1,num2str(newVals(jj)));
                    end
                end
                % Update Dragger radius
                dragInd = find(contains(inText(rb1:end),'<DraggerSize>'),1,'first')+rb1-1;
                dragSize = double(extractBetween(string(inText{dragInd}),'>','</'));
                if dragSize ~= -1
                    inText{dragInd} = replaceBetween(inText{dragInd},'>','</',num2str(factor*dragSize));
                end
                % Update Attachment Radius
                if contains(inText{typeInd},'Attachment')
                    radiusInd = find(contains(inText(rb1:end),'<Radius>'),1,'first')+rb1-1;
                    radiusSize = double(extractBetween(string(inText{radiusInd}),'>','</'));
                    inText{radiusInd} = replaceBetween(inText{radiusInd},'>','</',num2str(factor*radiusSize));
                end
            else % If muscle
                % Change the muscle parameters
                ksInd = find(contains(inText(rb1:end),'<Kse>'),1,'first')+rb1-1;
                kpInd = find(contains(inText(rb1:end),'<Kpe>'),1,'first')+rb1-1;
                bInd = kpInd+1; % This is easier than searching for the second B Value (since the ST curve uses a B Value, too)                
                muscParamInds = [ksInd,kpInd,bInd];
                if nargin < 3 % If less than three inputs
                    % No forced VE params, just act normally
                    for jj = 1:length(muscParamInds)
                        oldVal = double(extractBetween(string(inText{muscParamInds(jj)}),'>','</'));
                        oldVals(jj) = oldVal;
                        newVEs(jj) = factor*oldVal;
                        inText{muscParamInds(jj)} = replaceBetween(inText{muscParamInds(jj)},'>','</',num2str(factor*oldVal));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         inText{muscParamInds(jj)} = replaceBetween(inText{muscParamInds(jj)},'>','</',num2str(1*oldVal));
%                         newVEs(jj) = 1*oldVal;
                    end
                else
                    % Force VE params given, used to test scaling
                    inText{muscParamInds(1)} = replaceBetween(inText{muscParamInds(1)},'>','</',num2str(inVals(1)));
                    inText{muscParamInds(2)} = replaceBetween(inText{muscParamInds(2)},'>','</',num2str(inVals(2)));
                    inText{muscParamInds(3)} = replaceBetween(inText{muscParamInds(3)},'>','</',num2str(inVals(3)));
                end
                % For an individual Fmax scale. Can eliminate next two lines and add fmInd to muscParamInds to revert this action
                    fmInd = find(contains(inText(rb1:end),'<MaximumTension>'),1,'first')+rb1-1;
                    oldFM = double(extractBetween(string(inText{fmInd}),'>','</'));
                    if factor < 1
                        newFM = factor^2*oldFM;
                    else
                        newFM = factor^3*oldFM;
                    end
                    inText{fmInd} = replaceBetween(inText{fmInd},'>','</',num2str(newFM));
                    %inText{fmInd} = replaceBetween(inText{fmInd},'>','</',num2str(100000));
                % Change LT params
                    ltInd = find(contains(inText(rb1:end),'<LengthTension>'),1,'first')+rb1-1;
                    lrInd = find(contains(inText(ltInd:end),'RestingLength'),1,'first')+ltInd-1; % Resting Length
                    lwInd = find(contains(inText(ltInd:end),'Lwidth'),1,'first')+ltInd-1; % Lwidth
                    llInd = find(contains(inText(ltInd:end),'LowerLimit'),1,'first')+ltInd-1; % Lower Limit
                    ulInd = find(contains(inText(ltInd:end),'UpperLimit'),1,'first')+ltInd-1; % Upper Limit
                    ltIndVec = [lrInd, lwInd, llInd, ulInd];
                    for jj = 1:length(ltIndVec)
                        oldVal = double(extractBetween(string(inText{ltIndVec(jj)}),'>','</'));
                        inText{ltIndVec(jj)} = replaceBetween(inText{ltIndVec(jj)},'>','</',num2str(factor*oldVal));
                    end
                % Change ST params
                    newSTmax = (1+newVEs(2)/newVEs(1))*newFM;
                    newYoff = -.01*newSTmax;
                    stInd = find(contains(inText(rb1:end),'<StimulusTension>'),1,'first')+rb1-1;
                    stMaxInd = find(contains(inText(stInd:end),'<B>'),1,'first')+stInd-1; % STmax
                    yOffInd = find(contains(inText(stInd:end),'<D>'),1,'first')+stInd-1; % Lwidth
                    inText{stMaxInd} = ['<B>',num2str(newSTmax),'</B>'];
                    inText{stMaxInd-2} = ['<UpperOutput>',num2str(newSTmax),'</UpperOutput>'];
                    inText{yOffInd} = ['<D>',num2str(newYoff),'</D>'];
            end
        end
        jIndVec = find(contains(inText,'<Joint>'));
        % Update JointSize
        for ii = 1:length(jIndVec)
            j1 = jIndVec(ii);
            jSizeInd = find(contains(inText(j1:end),'<Size>'),1,'first')+j1-1;
            jSize = double(extractBetween(string(inText{jSizeInd}),'>','</'));
            inText{jSizeInd} = replaceBetween(inText{jSizeInd},'>','</',num2str(factor*jSize));
        end
    end
    
    % Move floor down
        scaledYVal = getYVal(inText,is_proj,'<ID>1adc8754-4d4b-4b4f-8085-53b73bf2a758</ID>');
        projObj = FullLeg(inPath,[],[]);
        attachCell = projObj.musc_obj{19}.pos_attachments(6,:);
        attachCell{1} = [0;scaledYVal;0];
        toeYVal = projObj.att_pos_on_demand([0;0;0],attachCell);
        toeYVal = toeYVal(2);
        inText = setYVal(inText,is_proj,'67e38fba-f56d-4448-b9a9-0b5cc7c62789',toeYVal+factor*toeYVal);
    
    % Write new project file
        % define new project name
        [path2,inName,inExt] = fileparts(inPath);
        scaledPathName = [path2,'\',inName,'_scaled',inExt];
        if is_proj
            inText{2} = ['<ProjectName>',inName,'_scaled</ProjectName>'];
            inText{3} = ['<SimulationFile>',inName,'_scaled_Standalone.asim</SimulationFile>'];
        end
        % write to new name
        fileID = fopen(scaledPathName,'w');
        fprintf(fileID,'%s\n',inText{:});
        fclose(fileID);
end

function yVal = getYVal(inText,is_proj,objID)
    % Input: inText: cell array of project text
    % Input: objID: char: ID of object pulling Y position from
    objInd = find(contains(inText,objID));
    if is_proj
        yValLine = inText{find(contains(inText(objInd:end),'<Y Value="'),1,'first')+objInd-1};
        yVal = extractBetween(yValLine,'"','"');
        yVal = str2double(yVal{1});
    else
        yValLine = inText{find(contains(inText(objInd:end),'<Position'),1,'first')+objInd-1};
        yVal = double(extractBetween(string(yValLine),'y="','"'));
    end
end

function inText = setYVal(inText,is_proj,objID,yVal)
    % Input: inText: cell array of project text
    % Input: objID: char: ID of object to set Y position of
    % yVal: double: y value to set to
    objInd = find(contains(inText,objID));
    if is_proj
        yValInd = find(contains(inText(objInd:end),'<Y Value="'),1,'first')+objInd-1;
        yTemp = @(inVal) ['<Y Value="',num2str(inVal),'" Scale="None" Actual="',num2str(inVal),'"/>'];
        inText{yValInd} = yTemp(yVal);
        % Change local matrices
        temp = find(contains(inText(objInd:end),'LocalMatrix'),2,'first')+objInd-1;
        lm1Ind = temp(1);
        lm2Ind = temp(2);
        lmInds = [lm1Ind, lm2Ind];
        for ii = 1:length(lmInds)
            oldMat = sscanf(char(extractBetween(string(inText{lmInds(ii)}),'>','</')), '%g,')';
            oldMat(14) = 10*yVal;
            newMat = [];
            for jj = 1:length(oldMat)
                if oldMat(jj) == 0
                    newMat = [newMat,'0,'];
                elseif oldMat(jj) == 1
                    newMat = [newMat,'1,'];
                elseif oldMat(jj) == -1
                    newMat = [newMat,'-1,'];
                else
                    newMat = [newMat,[sprintf('%.6f',oldMat(jj)),',']];
                end
            end
            newMat = newMat(1:end-1);
            inText{lmInds(ii)} = replaceBetween(inText{lmInds(ii)},'>','</',newMat);
        end
    else
        yValInd = find(contains(inText(objInd:end),'<Position'),1,'first')+objInd-1;
        inText{yValInd} = replaceBetween(inText{yValInd},'y="','"',num2str(yVal));
    end
end

function newLine = replacePositionInLine(oldLine,newVal)
    temp1 = replaceBetween(oldLine,'Value="','"',num2str(newVal));
    temp2 = replaceBetween(temp1,'Actual="','"',num2str(newVal));
    newLine = replaceBetween(temp2,'Scale="','"','None');
end

function inText = changeWalkingWaveforms(inText,factor,is_proj)
    load([pwd,'\Data\walking_waveforms.mat'],'walking_waveforms')
    owp = walking_waveforms.NW_walking_motion; % original walking pattern
    ncp = .5519*factor.^.32; % new cycle period, relationship from Hooper 2009 for trot, .5519 was cycle period for rat walking found from NW data
    physT = str2double(extractBetween(string(inText{13}),'>','<')); % physics timestep of the system
    ntv = (0:physT:3*ncp)'; % new time vector
    nwp = interp1(1:length(owp),owp,linspace(1,length(owp),length(ntv)));
    % Fill in joint info
    stimNames = {'Walking_Hip','Walking_Knee','Walking_xAnkle'};
    for ii = 1:3
        [~, ~, ~, equations, timeBnds] = sumsinesFit(ntv, nwp(:,ii),8);
        stimInd = find(contains(inText,stimNames(ii)));
        if ~is_proj
            % Update stimulus
                inText{stimInd+2} = '<Enabled>True</Enabled>';
                inText{stimInd+8} = ['<EndTime>',num2str(timeBnds(2)),'</EndTime>'];
                inText{stimInd+10} = ['<Equation>',equations{2},'</Equation>'];
        else
        end
    end
        % Update sim end time
            inText{9} = ['<SimEndTime>',num2str(timeBnds(2)+.1),'</SimEndTime>'];
        % Update Datatool end times
        jmInd = find(contains(inText,'<Name>JointMotion'));
        inText{jmInd+5} = ['<EndTime>',num2str(timeBnds(2)),'</EndTime>'];
        ptInd = find(contains(inText,'<Name>PassiveTension'));
        inText{ptInd+5} = ['<EndTime>',num2str(timeBnds(2)),'</EndTime>'];
end