%[obj,sim_file,joints,bodies,joint_limits,joint_profile,sdata,passive_tension] = design_synergy("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_forceSensor_Standalone.asim");
simPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_forceSensor_auto.asim';
obj = design_synergy(simPath);

%% Input the muscle stimulus for the leg position trials
muscle_info = scrapeFileForMuscleInfo(simPath);

legPositions =  [-0.952,-0.169,-0.473;...
                -0.874,-0.293986,-0.585245;...
                -0.685765,-0.387926,-0.752462;...
                -0.395308,-0.385050,-0.855123;...
                -0.0815863,-0.277706,-0.89230;...
                0.127572,-0.143481,-0.724004;...
                0.132435,-0.228251,-0.452663;...
                -0.0692101,-0.517073,-0.515119;...
                -0.380566,-0.700281,-0.891764;...
                -0.673678,-0.683119,-1.12102;...
                -0.868397,-0.532870,-1.02805;...
                -0.958282,-0.294953,-0.747424;...
                -0.971552,-0.161160,-0.510542];
    
legPos2run = 5:7;%1:size(legPositions,1);
musc2run = 1:size(muscle_info,1);
forceData = cell(max(musc2run),max(legPos2run));
runCounter = 1;

for uu = musc2run
    stimMuscle = muscle_info{uu,2};

    simText = importdata(simPath);

    % Find the muscle's connected neuron name
    temp = find(contains(simText,['ad-',stimMuscle]),1,'first')+6;
    sConnectedNeuron = extractBetween(simText{temp},'>','<');

    % Modify the muscle stimulus object
    sInd = find(contains(simText,'692c0e30-fb6a-470b-94e6-2940eb4ee869'),1,'first');
    sName = find(contains(simText(sInd:end),'<Name>'),1,'first')+sInd-1;
    simText{sName} = replaceBetween(simText{sName},'>','</',[stimMuscle(4:end),'Stim']);
    sNeur = find(contains(simText(sInd:end),'<TargetNodeID>'),1,'first')+sInd-1;
    simText{sNeur} = replaceBetween(simText{sNeur},'>','</',sConnectedNeuron{1});
    
    jointStims = {'ConstantHip';'ConstantKnee';'ConstantxAnkle'};

    for kk = legPos2run
        position = legPositions(kk,:);
        for joint = 1:3
            jInd = find(contains(simText,jointStims{joint}),1,'first');
            jPosInd = find(contains(simText(jInd:end),'Equation'),1,'first')+jInd-1;
            simText{jPosInd} = ['<Equation>',num2str(position(joint)),'</Equation>'];
        end

        fid = fopen(simPath,'w');
        fprintf(fid,'%s\n',simText{:});
        fclose(fid);

        sdata = processSimData(simPath);

        initForcePos_rel = [0;-38.242;-9.6]./1000;
        anchorPos = [50.061;-139.07;18.55]./1000;

        temp = find(contains({sdata.name},'KeyMNs'));
        [~,activeNeur] = max(max(sdata(temp).data(1:end-10,:)));
        stimWave = sdata(temp).data(1:end-10,activeNeur);
        activeInds = find(stimWave>-.06);

        temp = find(contains({sdata.name},'JointMotion'));
        jointMotionAKH = sdata(temp).data;
        jointMotion = jointMotionAKH;
        jointMotion(:,1) = jointMotionAKH(:,3);
        jointMotion(:,3) = jointMotionAKH(:,1);
        
        initLegPos = jointMotion(1866,:);
        initForcePos = obj.att_pos_on_demand(initLegPos,[{initForcePos_rel},{''},{3},{''}]);
        temp1 = find(contains({sdata(:).name},'ForceSensorInfo'));
        temp2 = find(contains(sdata(temp1).data_cols,'Ib'));
        temp3 = find(contains(sdata(temp1).data_cols,'ml'));
        fsData = sdata(temp1).data(:,temp2);
        fsDataL = sdata(temp1).data(:,temp3);
        fsBase = mean(fsData(12000:16000));
        clear temp* fX fY fZ

        counter = 1;
        for ii = 1:length(jointMotion)%activeInds(1):activeInds(end)%floor(linspace(1,18529,50))
            movedForcePos = obj.att_pos_on_demand(jointMotion(ii,:),[{initForcePos_rel},{''},{3},{''}]);

            fVec = movedForcePos-initForcePos;
            if sum(fVec) == 0
                fVecnorm2(counter,:) = [0,0,0];
            else
                fVecnorm2(counter,:) = (fVec./norm(fVec))';
            end

            fAng_all(counter,1) = atan2(fVecnorm2(counter,2),-fVecnorm2(counter,1))*180/pi;

            forcePosLog(counter,:) = movedForcePos';

            counter = counter + 1;
        end

%         forcePosRel = forcePosLog-mean(forcePosLog(17000:18000,:));
        forcePosRel = forcePosLog-mean(forcePosLog(200:400,:));
        %forcePosRel(:,1) = -1.*forcePosRel(:,1);
        sumForce = sqrt(sum(forcePosRel.^2,2));
        
        stimBnds = [find(stimWave>-.06,1,'first'),find(stimWave>-.06,1,'last')];
        maxStim = floor(sum(stimBnds)/2);
        fAng_max = mean(fAng_all(maxStim-400:maxStim+400));
        
        % Store data
        outData = struct();
        outData.name = stimMuscle;
        outData.force_rel = forcePosRel;
        outData.force = forcePosLog;
        outData.force_sum = sumForce;
        outData.fAng_all = fAng_all;
        outData.fAng_max = fAng_max;
        outData.fVec_norm = fVecnorm2;
        outData.jointMotion = jointMotion;
        
        forceData{uu,kk} = outData;
        disp([stimMuscle,' (',num2str(uu),'/',num2str(max(musc2run)),') for leg position (',num2str(kk-legPos2run(1)+1),'/',num2str(length(legPos2run)),') done. (',num2str(runCounter),'/',num2str(length(musc2run)*length(legPos2run)),')'])
        runCounter = runCounter + 1;
    end
    
end
%% Plotting Quivers
if 1
    [numMusc,numPos] = size(forceData);
    figure('Position',[962,45,958,950])
    for ii = musc2run
        subObjs(ii) = subplot(ceil(numMusc/4),4,ii);
        for jj = legPos2run
            fVecData = forceData{ii,jj}.fVec_norm;
            [~,maxfVecInd] = max(sum(fVecData,2));
            fvec = fVecData(maxfVecInd,:);
            quiver(0,0,fvec(1),fvec(2))
            fbnds(ii,:) = fvec(1:2);
            hold on
            axis('equal')
            title(muscle_info{ii,2}(4:end))
        end
    end
    fbnds(end+1,:) = [-.1,-.1];
    xlim(subObjs,[min(fbnds(:,1)), max(fbnds(:,1))]);
    ylim(subObjs,[min(fbnds(:,2)), max(fbnds(:,2))]);
end
%% Plotting Quivers for All Positions
if 0
    mnum = 4;
    figure('Position',[962,45,958,950])
    neutralLeg = [0,0,0];
    % Plot the neutral leg, muscle, and force sensor positions
        movedForcePos = obj.att_pos_on_demand(neutralLeg,[{initForcePos_rel},{''},{3},{''}]);
        forcePosLog(counter,:) = movedForcePos';

        pospel = obj.organism_position*1000;
        posfem = obj.att_pos_on_demand(neutralLeg,[{[0;0;0]},{''},{2},{''}])'*1000;
        postib = obj.att_pos_on_demand(neutralLeg,[{[0;0;0]},{''},{3},{''}])'*1000;
        posfot = obj.att_pos_on_demand(neutralLeg,[{[0;0;0]},{''},{4},{''}])'*1000;
        postoe = obj.att_pos_on_demand(neutralLeg,[{[20.452;-6.26;-1.398]./1000},{''},{4},{''}])'*1000;
        limbLine = [pospel;posfem;postib;posfot;postoe];
        plot3(limbLine(:,1),limbLine(:,2),limbLine(:,3),'LineWidth',3)
        hold on
        plot3(movedForcePos(1)*1000,movedForcePos(2)*1000,movedForcePos(3)*1000,'r.','MarkerSize',20)
        muscMat = zeros(size(obj.musc_obj{mnum}.pos_attachments,1),3);
        for jj = 1:size(obj.musc_obj{mnum}.pos_attachments,1)
            muscMat(jj,:) = obj.att_pos_on_demand(neutralLeg,obj.musc_obj{mnum}.pos_attachments(jj,:))*1000;
        end
        plot3(muscMat(:,1),muscMat(:,2),muscMat(:,3),'m','LineWidth',3)
        grid on
        axis('equal')
        xlabel('X');ylabel('Y');zlabel('Z');
        xlim([-50 40]);ylim([-180 -90]);zlim([0 25]);
        view([0 90])
        title([forceData{mnum,legPos2run(1)}.name(4:end)],'FontSize',20)
    % Go through the different positions and plot the quivers
        for ii = 1:length(legPos2run)
            legPos = legPositions(legPos2run(ii),:);
            fPosTemp = obj.att_pos_on_demand(legPos,[{initForcePos_rel},{''},{3},{''}]).*1000; 
            plot3(fPosTemp(1),fPosTemp(2),fPosTemp(3),'k.','MarkerSize',20)
            hold on
            fVecData = forceData{mnum,legPos2run(ii)}.fVec_norm;
            % WHICH QUIVER DO I PLOT??
            [~,maxfVecInd] = max(sum(fVecData(500:end,:),2));
            %fVDtemp = fVecData(3000:6000,:);
            %[~,maxfVecInd] = max(diff(sum(fVDtemp.^2,2)));
            %[~,maxfVecInd] = max(sum(fVDtemp.^2,2));
            maxfVecInd = maxfVecInd+500;
            fvec = fVecData(maxfVecInd,:);
            quiver(fPosTemp(1),fPosTemp(2),10*fvec(1),10*fvec(2),'r','LineWidth',1.5,'MaxHeadSize',.5)
        end
end
%% Plotting Motion Over Time
if 0
    dataCell = forceData{6,8};
    close all
    figure('Position',[962,45,958,950])
    counter = 1;
    initForcePos = mean(dataCell.force(200:3000,:));
    mnum = find(contains(muscle_info(:,2),dataCell.name));
    for ii = floor(linspace(activeInds(1)-1500,activeInds(end)+1500,100))%floor(linspace(1,18529,50))
        movedForcePos = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),[{initForcePos_rel},{''},{3},{''}]);
        forcePosLog(counter,:) = movedForcePos';

        %fVec = movedForcePos-anchorPos;
        fVec = movedForcePos-initForcePos;
        fVecnorm = dataCell.fVec_norm;
        fsData = 1000*dataCell.force_rel;

        fX(counter) = fsData(ii,1);
        fY(counter) = fsData(ii,2);
        fZ(counter) = fsData(ii,3);

        %fAng_all(counter) = atan2(fVecnorm(2),-fVecnorm(1))*180/pi;
        fAng_all = dataCell.fAng_all;

        pospel = obj.organism_position*1000;
        posfem = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),[{[0;0;0]},{''},{2},{''}])'*1000;
        postib = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),[{[0;0;0]},{''},{3},{''}])'*1000;
        posfot = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),[{[0;0;0]},{''},{4},{''}])'*1000;
        postoe = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),[{[20.452;-6.26;-1.398]./1000},{''},{4},{''}])'*1000;
        limbLine = [pospel;posfem;postib;posfot;postoe];
        plot3(limbLine(:,1),limbLine(:,2),limbLine(:,3))
        hold on
        plot3(movedForcePos(1)*1000,movedForcePos(2)*1000,movedForcePos(3)*1000,'r.','MarkerSize',20)
        switch 2
            case 1 % For plotting individual x,y,z components of force vector
                xMat = [movedForcePos'*1000;movedForcePos'*1000+[fX(counter),0,0]*15];
                yMat = [movedForcePos'*1000;movedForcePos'*1000+[0,fY(counter),0]*15];
                zMat = [movedForcePos'*1000;movedForcePos'*1000+[0,0,fZ(counter)]*15];
                plot3(xMat(:,1),xMat(:,2),xMat(:,3),'k','LineWidth',2)
                plot3(yMat(:,1),yMat(:,2),yMat(:,3),'b','LineWidth',2)
                plot3(zMat(:,1),zMat(:,2),zMat(:,3),'g','LineWidth',2)
            case 2 % For plotting the combined force vector
                fMat = [movedForcePos'*1000;movedForcePos'*1000+[fX(counter),fY(counter),fZ(counter)]*10];
                plot3(fMat(:,1),fMat(:,2),fMat(:,3))
            case 3 % For xy and yz planes of force vector
                xyMat = [movedForcePos'*1000;movedForcePos'*1000+[fX(ii),fY(ii),0]*10];
                yzMat = [movedForcePos'*1000;movedForcePos'*1000+[0,fY(ii),fZ(ii)]*10];
                plot3(xyMat(:,1),xyMat(:,2),xyMat(:,3),'k','LineWidth',3)
        end
        muscMat = zeros(size(obj.musc_obj{mnum}.pos_attachments,1),3);
        for jj = 1:size(obj.musc_obj{mnum}.pos_attachments,1)
            %muscMat(jj,:) = obj.musc_obj{mnum}.pos_attachments{jj,4}(ii,:)*1000;
            muscMat(jj,:) = obj.att_pos_on_demand(dataCell.jointMotion(ii,:),obj.musc_obj{mnum}.pos_attachments(jj,:))*1000;
        end
        plot3(muscMat(:,1),muscMat(:,2),muscMat(:,3),'m')
        plot3(initForcePos(1)*1000,initForcePos(2)*1000,initForcePos(3)*1000,'k.','MarkerSize',10)
        %plot3(yzMat(:,1),yzMat(:,2),yzMat(:,3),'b')
        grid on
        axis('equal')
        xlabel('X');ylabel('Y');zlabel('Z');
        xlim([-50 40]);ylim([-180 -90]);zlim([0 25]);
        view([0 90])
        title(num2str(obj.theta_motion_time(ii,1)),'FontSize',20)
        pause(.1)
        hold off
        counter = counter + 1;
    end
end
%% Plotting section force waveforms
if 0
    figure('Position',[680,35,781,943]);
    subplot(3,1,1)
    plot(forcePosRel)
            xlabel('Time (ms)')
            ylabel('Force')
            title('BFP Force XYZ')
            legend({'X';'Y';'Z'})
            %xlim([0 10*length(mkdat(:,33))])
    subplot(3,1,2)
    plot(sumForce)
            xlabel('Time (ms)')
            ylabel('Force')
            title('BFP Sum Force')
            %xlim([0 10*length(mkdat(:,33))])
    subplot(3,1,3)
    plot(fAng_all)
            xlabel('Time (ms)')
            ylabel('F Angle')
            xlim([0 length(fAng_all)])
end