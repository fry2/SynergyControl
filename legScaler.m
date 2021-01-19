function legScaler(normal_project_path)
    % For an input project path, scale the leg up by a certain factor
    % Input: normal_project_path (char): full path to the project to be scaled
    factor = 2;
    
    % Import project text
    if ~ischar(normal_project_path)
        normal_project_path = char(normal_project_path);
    end
    projtext = importdata(normal_project_path);
    
    % Find organism position
        tempind = find(contains(projtext,'<Organism>'))+17;
        counter = 1; orgpos = zeros(1,3);
        for ii = tempind+1:tempind+3
            orgpos(counter) = double(extractBetween(string(projtext{ii}),'Actual="','"'));
            counter = counter + 1;
        end
    
    % Update all scaling factors in the project text
    scaleIndVec = find(contains(projtext,'<Scale>'));
    for ii = 1:length(scaleIndVec)
        scaleInd = scaleIndVec(ii);
        for jj = scaleInd+1:scaleInd+3
            projtext{jj} = replacePositionInLine(projtext{jj},factor);
        end
    end
    
    % Update all local positions based on new scale
    rbIndVec = find(contains(projtext,'<RigidBody>'));
    for ii = 1:length(rbIndVec)
        rb1 = rbIndVec(ii);
        if ~contains(projtext{rb1+3},'<Type>LinearHillMuscle</Type>')
            lpInd = find(contains(projtext(rb1:end),'<LocalPosition>'),1,'first')+rb1-1;
            counter = 1;
            for jj = lpInd+1:lpInd+3
                oldVal = double(extractBetween(string(projtext{jj}),'Actual="','"'));
                projtext{jj} = replacePositionInLine(projtext{jj},factor*oldVal);
                newlp(counter) = factor*oldVal;
                counter = counter + 1;
            end
            lpMatInd = find(contains(projtext(rb1:end),'<LocalMatrix'),1,'first')+rb1-1;
            oldMat = sscanf(char(extractBetween(string(projtext{lpMatInd}),'>','</')), '%g,')';
            if all(oldMat(13:15)==[0,0,0])
                newMat = [sprintf('%.0f,',oldMat(1:12)),sprintf('%.0f,' , newlp),sprintf('%.0f,',oldMat(16:end))];
            else
                newMat = [sprintf('%.0f,',oldMat(1:12)),sprintf('%.5f,' , newlp),sprintf('%.0f,',oldMat(16:end))];
            end
            newMat = newMat(1:end-1);
            projtext{lpMatInd} = replaceBetween(projtext{lpMatInd},'>','</',newMat);
        end
    end
    
    % Update bone densities
%     boneNames = {'Pelvis','LH_Femur','LH_Tibia','LH_Foot'};
%     for ii = 1:length(boneNames)
%         boneInd = find(contains(projtext,['<Name>',boneNames{ii},'</Name>']));
%         densInd = find(contains(projtext(boneInd:end),'<Density'),1,'first')+boneInd-1;
%         oldVal = double(extractBetween(string(projtext{densInd}),'Actual="','"'));
%         projtext{densInd} = replacePositionInLine(projtext{densInd},(1/factor^3)*oldVal);
%     end
    
    % Update bone masses
    boneNames = {'Pelvis','LH_Femur','LH_Tibia','LH_Foot'};
    for ii = 1:length(boneNames)
        boneInd = find(contains(projtext,['<Name>',boneNames{ii},'</Name>']));
        massInd = find(contains(projtext(boneInd:end),'<Mass'),1,'first')+boneInd-1;
        mass = double(extractBetween(string(projtext{massInd}),'Actual="','"'));
        projtext{massInd} = replacePositionInLine(projtext{massInd},factor^3*mass);
    end
    
    % Move floor down
        scaledYVal = getYVal(projtext,'<ID>1adc8754-4d4b-4b4f-8085-53b73bf2a758</ID>');
        projObj = FullLeg(normal_project_path,[],[]);
        attachCell = projObj.musc_obj{19}.pos_attachments(6,:);
        attachCell{1} = [0;scaledYVal;0];
        toeYVal = projObj.att_pos_on_demand([0;0;0],attachCell);
        toeYVal = toeYVal(2);
        projtext = setYVal(projtext,'67e38fba-f56d-4448-b9a9-0b5cc7c62789',toeYVal+factor*toeYVal);
    
    % Write new project file
        % define new project name
        [projPath,projName,projExt] = fileparts(normal_project_path);
        scaledProjectName = [projPath,'\',projName,'_scaled',projExt];
        projtext{2} = ['<ProjectName>',projName,'_scaled</ProjectName>'];
        projtext{3} = ['<SimulationFile>',projName,'_scaled_Standalone.asim</SimulationFile>'];
        % write to new name
        fileID = fopen(scaledProjectName,'w');
        fprintf(fileID,'%s\n',projtext{:});
        fclose(fileID);
end

function yVal = getYVal(projtext,objID)
    % Input: projtext: cell array of project text
    % Input: objID: char: ID of object pulling Y position from
    objInd = find(contains(projtext,objID));
    yValLine = projtext{find(contains(projtext(objInd:end),'<Y Value="'),1,'first')+objInd-1};
    yVal = extractBetween(yValLine,'"','"');
    yVal = str2double(yVal{1});
end

function projtext = setYVal(projtext,objID,yVal)
    % Input: projtext: cell array of project text
    % Input: objID: char: ID of object to set Y position of
    % yVal: double: y value to set to
    objInd = find(contains(projtext,objID));
    yValInd = find(contains(projtext(objInd:end),'<Y Value="'),1,'first')+objInd-1;
    %yValLine = projtext{yValInd};
    %yLineOut = replaceBetween(yValLine,'Value="','"',num2str(yVal));
    yTemp = @(inVal) ['<Y Value="',num2str(inVal),'" Scale="None" Actual="',num2str(inVal),'"/>'];
    projtext{yValInd} = yTemp(yVal);
end

function newLine = replacePositionInLine(oldLine,newVal)
%     scaleVec = [contains(oldLine,'None'),contains(oldLine,'pico'),contains(oldLine,'nano'),contains(oldLine,'milli'),...
%     contains(oldLine,'centi'),contains(oldLine,'Kilo'),contains(oldLine,'Giga')];
%     scaleVals = [1,1e-12,1e-9,1e-3,1e-2,1e3,1e9];
    temp1 = replaceBetween(oldLine,'Value="','"',num2str(newVal));
    temp2 = replaceBetween(temp1,'Actual="','"',num2str(newVal));
    newLine = replaceBetween(temp2,'Scale="','"','None');
end