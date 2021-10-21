inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_walking.asim";
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_Size_Rat_Walk_Rat.asim";
%inSimPath = "C:\Users\fry16\OneDrive\Documents\JointDampingOpt\InjectedProject\JointDampingOpt_injected_Standalone.asim";

    if contains(inSimPath,'Walk_Cat')
        beg = 385; mid = 724; ennd = 1056; % is cat
    elseif contains(inSimPath,'Walk_Horse')
        beg = 895; mid = 1479; ennd = 2300; % is horse
    else
        beg = 654; mid = 1034; ennd = 1676;  % is else (rat)
    end 
    oneHundredizer = @(inMat,trans) [interp1(1:length(inMat(1:trans,:)),inMat(1:trans,:),linspace(1,length(inMat(1:trans,:)),37)),...
                                     interp1(1:length(inMat(trans:end,:)),inMat(trans:end,:),linspace(1,length(inMat(trans:end,:)),63))];

lengthVals = exp(linspace(log(.05),log(25),24)); 
if ~any(lengthVals==1)
    ind = find((1-lengthVals)>0,1,'last');
    lengthVals = [lengthVals(1:ind),1,lengthVals(ind+1:end)]; 
end
%lengthVals = [.05,.36,.75,1,1.3,12]; % SWING CROSSOVERS
%lengthVals = [.05,.215,.53,.615,1,12]; % STANCE CROSSOVERS
%lengthVals = [.05,.25,.5,1,2,3,7,9,12]; 
%lengthVals = [.375,.9,1,2.63,6.46,12]; % STANCE ANIMALS
%lengthVals = [2.63,6.46,8.66]; % STANCE ANIMALS
%lengthVals = [.05,.2,.5,1,3,12];
numLens = length(lengthVals);
% Mean rate change of moment arms
moMean = @(mo) mean(abs(diff(mo(sum(mo(1:38,:),2)~=0,100:end-100)')./.54e-3),'all');
count = 1;

% run an initial test to check for sim data sizing
obj = design_synergy(inSimPath);
dat_len = length(obj.theta_motion); 

% Pre-allocate for speed
% elasK = zeros(18529,length(lengthVals),3); gravK = elasK; elasK_approx = elasK; gravK_approx = elasK; bRot = elasK; bRot_approx = elasK;
% jmDoc = []; moDoc = zeros(length(lengthVals),3); zeta = jmDoc; 
force_mn = cell(1,length(lengthVals)); force_tot = force_mn;
body_masses = zeros(4,length(lengthVals)); body_lengths = body_masses; body_densities = body_masses;

 exitflags = zeros(dat_len,length(lengthVals)); objCell = cell(1,numLens); 
 torques_active = zeros(length(obj.theta_motion),3,numLens); 
 torques_active = objCell;
 
 torques_muscle = torques_active; torques_inertial = torques_active; jointmotion = torques_active; torques_load = torques_active;
 torques_grav = torques_active;

% Cycle through length scales for data
for ii = 1:length(lengthVals)
    tstart = tic;
    lengthScale = lengthVals(ii);
    outPath = legScaler(inSimPath,lengthScale);
        obj = design_synergy(outPath);
        [pks,locs] = findpeaks(obj.theta_motion(:,3)); 
        rangeMat(ii,:) = locs(2:4);
    % Gather data on force production during normal motion
        jointmotion{ii} = obj.theta_motion.*(180/pi)+[98.4373 102.226 116.2473];
        %results_cell = pedotti_optimization(obj);
        [torques_active{ii},torques_muscle{ii},torques_inertial{ii}] = compute_active_joint_torque(obj,0);
        torques_load{ii} = compute_load_torques(obj,0);
        torques_grav{ii} = compute_grav_torques(obj,0);
        %force_mn{ii} = results_cell{3,2};
        %force_tot{ii} = results_cell{2,2};
        %exitflags(:,ii) = results_cell{9,2};
        for bodyNum = 1:4
            body = obj.body_obj{bodyNum};
            body_masses(bodyNum,ii) = body.mass;
            body_densities(bodyNum,ii) = body.density;
            if bodyNum > 1
                body_lengths(bodyNum-1,ii) = body.length;
            end
        end
        objCell{ii} = obj;
    % Display and increment
    telapsed = toc(tstart);
    if telapsed < 60
        timeStr = ['(',num2str(telapsed),' sec)'];
    else
        timeStr = ['(',num2str(telapsed/60),' min)'];
    end
    disp(['Completed ',num2str(count),' out of ',num2str(length(lengthVals)),'. ',timeStr])
    count = count +1;
end
% Zoning Sorter
muscZones = zoning_sorter(inSimPath,6); mDataZone = cell(6,1); am_all = cell(1,length(lengthVals));torques_act = cell(1,3); torques_tot = cell(1,3);
for jj = 1:length(lengthVals)
    mTemp = cell(6,1);
    %[Am_musc,V_musc,Am_musc_raw] = Am_generator(objCell{jj},force_tot{jj});
    %am_all{jj} = Am_musc';
    for ii = 1:38
        muscName = muscZones{ii,1}; muscZone = muscZones{ii,2};
        if isempty(mTemp{muscZone,1})
            colNum = 1;
        else
            colNum = size(mTemp{muscZone,1},2)+1;
        end
%         mTemp{muscZone,1}(:,colNum) = force_mn{jj}(:,ii);
        %mTemp{muscZone,1}(:,colNum) = Am_musc(ii,:)';
        mTemp{muscZone,2}(:,colNum) = {muscName};
    end
    mDataZone{jj} = mTemp;  
    rm = find_relevant_muscles(objCell{1});
    for ii = 1:3
        %[moment_output] = compute_joint_moment_arms(objCell{jj},ii,1);
        %moment_output = moment_output(rm{ii},:)'./1000;
        %torque_temp = sum(moment_output.*force_mn{jj}(:,rm{ii}),2);
        %torques_act{ii}(:,jj) = normalizeTorque(torque_temp(beg:ennd),mid-beg);
        beg = rangeMat(jj,1); mid = rangeMat(jj,2); ennd = rangeMat(jj,3);
        torques_act{ii}(:,jj) = -normalizeTorque(squeeze(torques_active{jj}(beg:ennd,ii)),mid-beg,0);
        torques_full{ii}(:,jj) = -normalizeTorque(squeeze(torques_active{jj}(beg:ennd,ii)),mid-beg,1);
    end
end
%% Plot individual zones across the length scale
aa = figure('Position',[962,2,958,994]);
    zoneNum = 3;  %numMusc = size(mDataZone{lengthNum}{zoneNum},2); 
    stepInds = 654:1676; %[654,1034,1676,2057]
    timeVec = obj.theta_motion_time(stepInds)-obj.theta_motion_time(stepInds(1));
    zoneNames = {'HipAdductors';'HipFlexors';'KneeExtensors';'KneeFlexors';'AnklePlantarflexors';'AnkleDorsiflexors'};
    savePath = ['G:\My Drive\Rat\MeetingFiles\Meeting_20210628\',zoneNames{zoneNum},'\'];
    if exist(savePath, 'dir') ~= 7
       mkdir(savePath)
    end
for lengthNum = 1:numLens
    subplot(2,1,1)
%         plot(timeVec,jointmotion(stepInds,:,lengthNum),'LineWidth',3); grid on
%         ylabel('Joint angle (deg)');
        plot(timeVec,torques_active(stepInds,:,lengthNum),'LineWidth',3); grid on
        ylabel('Joint Torque (Nm)','FontSize',14);
        legend({'Hip';'Knee';'Ankle'},'Location','southwest','FontSize',10); xlim([min(timeVec) max(timeVec)])
        title([zoneNames{zoneNum},'. Length scale: ',num2str(lengthVals(lengthNum),3),'x rat. (',num2str(.2*lengthVals(lengthNum),3),' m)'],'FontSize',16)
        xlabel('Time (s)','FontSize',14);
    subplot(2,1,2)
        data = mDataZone{lengthNum}{zoneNum,1}(stepInds,:);
        plot(timeVec,data,'LineWidth',3); grid on
        legend(mDataZone{lengthNum}{zoneNum,2},'Interpreter','none','Location','northeast','FontSize',10)
        ylabel('Muscle Activation (N)','FontSize',14); xlabel('Time (s)','FontSize',14); 
        ylim([min(data,[],'all') max(data,[],'all')]); xlim([min(timeVec) max(timeVec)])
        xline(timeVec(1034-654),'HandleVisibility','off')
     drawnow
%        %saveas(aa,[savePath,'\Swing\',zoneNames{zoneNum},'_',num2str(lengthNum),'.png'])
%       saveas(aa,[savePath,zoneNames{zoneNum},'_',num2str(lengthNum),'.png'])
%       % Capture the plot as an image
%       frame = getframe(aa); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if lengthNum == 1 
%           imwrite(imind,cm,[savePath,zoneNames{zoneNum},'.gif'],'gif', 'Loopcount',inf,'DelayTime',2); 
%       else 
%           imwrite(imind,cm,[savePath,zoneNames{zoneNum},'.gif'],'gif','WriteMode','append','DelayTime',2); 
%       end 
end
%% Plot all total force profile across length scale GIF
aa = figure('Position',[962,2,958,994]); 
savePath = 'G:\My Drive\Rat\MeetingFiles\Meeting_20210628\TotalForce\';
    if exist(savePath, 'dir') ~= 7
       mkdir(savePath)
    end
for lengthNum = 1:length(lengthVals)
    numMusc = size(mDataZone{lengthNum}{zoneNum},2); 
    stepInds = 1:length(force_mn{1});stepInds = beg:ennd;timeVec = obj.theta_motion_time(stepInds)-obj.theta_motion_time(1);
    zoneNames = {'Hip adductors';'Hip flexors';'Knee extensors';'Knee flexors';'Ankle plantarflexors';'Ankle dorsiflexors'};

    subplot(2,1,1)
        plot(timeVec,torques_active(stepInds,:,lengthNum),'LineWidth',3); grid on
        ylabel('Joint Torque (Nm)','FontSize',14);
        legend({'Hip';'Knee';'Ankle'},'Location','southwest','FontSize',10); 
        xlim([min(timeVec) max(timeVec)])
        title([zoneNames{zoneNum},'. Length scale: ',num2str(lengthVals(lengthNum),3),'x rat. (',num2str(.2*lengthVals(lengthNum),3),' m)'],'FontSize',16)
        xlabel('Time (s)','FontSize',14);
    subplot(2,1,2)
        plot(timeVec,force_tot{lengthNum},'LineWidth',3); grid on
        %legend(mDataZone{lengthNum}{zoneNum,2},'Interpreter','none','Location','northwest','FontSize',10)
        ylabel('Muscle Force (N)'); xlabel('Time (s)'); xlim([min(timeVec) max(timeVec)]); ylim([0 max(force_tot{lengthNum},[],'all')]);
    drawnow
      saveas(aa,[savePath,'TotalForce_',num2str(lengthNum),'.png'])
      % Capture the plot as an image
      frame = getframe(aa); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if lengthNum == 1 
          imwrite(imind,cm,[savePath,'TotalForce','.gif'],'gif', 'Loopcount',inf,'DelayTime',2); 
      else 
          imwrite(imind,cm,[savePath,'TotalForce','.gif'],'gif','WriteMode','append','DelayTime',2); 
      end 
end
%% Plot individual muscle waveforms for the different scales
mNum = 12;
figure('Position',[962,2,958,994]); lengthVals2plot = [1,5,9,11];
normalized = 1; cm = cool(length(lengthVals)); r = linspace(0,1,length(lengthVals2plot))'; b = zeros(length(lengthVals2plot),1); 
g = linspace(1,0,length(lengthVals2plot))'; cm = [r,g,b];
stepInds = [654,1034,1676]; timeVec = obj.theta_motion_time(stepInds(1):stepInds(3))-obj.theta_motion_time(stepInds(1));
for ii = length(lengthVals2plot):-1:1
    lVal = lengthVals2plot(ii);
    if normalized
        plot(timeVec,am_all{lVal}(stepInds(1):stepInds(3),mNum)./max(am_all{lVal}(stepInds(1):stepInds(3),mNum)),'LineWidth',2,'Color',[cm(ii,:),1]);hold on
        ylabel('Muscle Activation'); xline(timeVec(stepInds(2)-stepInds(1)),'b','LineWidth',3,'HandleVisibility','off')
        title([obj.musc_obj{mNum}.muscle_name(4:end),' Activation v. Size Scale'],'FontSize',16)
    else
        semilogy(obj.theta_motion_time,force_tot{lVal}(:,mNum),'LineWidth',2);hold on
    end
end
legend(cellfun(@num2str,num2cell(fliplr(round(lengthVals(lengthVals2plot),2))),'UniformOutput',false))
xlabel('Time (s)'); xlim([min(timeVec) max(timeVec)])
%% Plot muscle group mean on semilogy plots
figure('Position',[962,2,958,994]);
cm = [0,1,0;1,0,0;0.9216,0.8235,0.2039;0,0,1;1,0,1;0,1,1];
zoneNames = {'HipAdductors';'HipFlexors';'KneeExtensors';'KneeFlexors';'AnklePlantarflexors';'AnkleDorsiflexors'};
for ii = 1:6
    for jj = 1:length(lengthVals)
        subplot(6,1,ii)
        meanData = mean(mDataZone{jj}{ii},2);
        meanData(meanData < 0) = 0;
        semilogy(obj.theta_motion_time,meanData,'Color',cm(ii,:),'Linewidth',2); hold on
        ylabel(zoneNames{ii}); xlim([0 max(obj.theta_motion_time)])
    end
end 
%% Plot normalized [TOTAL OR MUSCLE OR INERTIAL] joint torques across scale
figure('Position',[962,2,958,994]);
mode = 3;
switch mode
    case 1 % Total joint torque
        torque_data = torques_active;
        titleType = 'Total';
    case 2 % Passive muscle torque
        torque_data = torques_muscle;
        titleType = 'Passive';
    case 3 % Inertial torque
        torque_data = torques_inertial;
        titleType = 'Inertial';
end

stepInds = beg:ennd; cm = cool(length(lengthVals)); legendNames ={'Hip';'Knee';'Ankle'};

for jj = 1:3
    subplot(3,1,jj);
    tData = squeeze((torque_data(stepInds,jj,:) - min(torque_data(stepInds,jj,:))) ./ ( max(torque_data(stepInds,jj,:)) - min(torque_data(stepInds,jj,:)) ));
    tData = interp1(1:length(tData),tData,linspace(1,length(tData),100));
    for ii = 1:length(lengthVals)
        plot(tData(:,ii),'LineWidth',2,'Color',[cm(ii,:),1]); hold on;
    end
    ylabel('Normalized Joint Torque'); xlabel('Stride (%)'); xlim([0 100]); xline(37)
    colormap('cool'); cBar = colorbar; cBar.Ticks = linspace(0,1,5); cBar.TickLabels = num2cell(round(linspace(lengthVals(1),lengthVals(end),5),1));
    title([titleType,' ',legendNames{jj},' Torque'],'FontSize',16)
end
%% Plot Relative contributions of Muscle and Inertia across scale
    range = 4; jointNames = {'Hip';'Knee';'Ankle'}; corrIn = zeros(3,numLens); corrMs = corrIn; 
    switch range
        case 1
            range2plot = 654:780; % EARLY SWING
            phaseString = 'Early Swing';
        case 2
            range2plot = 780:906; % MID SWING
            phaseString = 'Mid Swing';
        case 3
            range2plot = 906:1034; % LATE SWING
            phaseString = 'Late Swing';
        case 4
            range2plot = 654:1034; % FULL SWING
            phaseString = 'Full Swing';
        case 5
            range2plot = 1034:1248; % EARLY STANCE
            phaseString = 'Early Stance';
        case 6
            range2plot = 1248:1462; % MID STANCE
            phaseString = 'Mid Stance';
        case 7
            range2plot = 1462:1676; % LATE STANCE
            phaseString = 'Late Stance';
        case 8
            range2plot = 1034:1676; % FULL STANCE
            phaseString = 'Full Stance';
    end
    figure('Position',[542,2,1378,994]);
    for joint = 1:3
        for len = 1:numLens
            range2plot = range_finder(range,len,rangeMat);
            corrIn(joint,len) = corr(torques_active{len}(range2plot,joint),torques_inertial{len}(range2plot,joint));
            corrMs(joint,len) = corr(torques_active{len}(range2plot,joint),torques_muscle{len}(range2plot,joint));
            corrGv(joint,len) = corr(torques_active{len}(range2plot,joint),torques_grav{len}(range2plot,joint));
            if range >=5
                corrLd(joint,len) = corr(torques_active{len}(range2plot,joint),torques_load{len}(range2plot,joint));
            end
        end
        [x0,y0] = find_crossover_point(corrMs,corrIn,lengthVals);
        subplot(3,1,joint) 
            semilogx([x0 x0],[y0-.5 y0+.5],'Color',[.75,.75,.75],'LineWidth',2); hold on
            semilogx(lengthVals,corrIn(joint,:),'r','LineWidth',2);
            semilogx(lengthVals,corrMs(joint,:),'b','LineWidth',2);
%             semilogx(lengthVals,corrGv(joint,:),'m','LineWidth',2);
%             if range >= 5
%                 semilogx(lengthVals,corrLd(joint,:),'g','LineWidth',2);
%             end
            ylim([-1.5,1.5]); grid on; ylabel('Relative Contribution','FontSize',16); xlabel('x Size Relative to Rat','FontSize',12);
            set(gca,'YTick',-1:.5:1)
            xlim([lengthVals(1),lengthVals(end)]);
            text(0.052,1.25,1,'Viscoelastic Forces','Color','b','FontSize',15);
            text(0.052,-1.25,1,'Inertial Forces','Color','r','FontSize',15);
        % Add animal labels
            xline(.375,'g-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
            xline(1,'c-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
            xline(2.63,'m-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
            xline(12,'y-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
            text(.173,-1.25,1,'Stick Insect \rightarrow','FontSize',15); 
            text(.716,-1.25,1,'Rat \rightarrow','FontSize',15); 
            text(1.89,-1.25,1,'Cat \rightarrow','FontSize',15); 
            text(7.66,-1.25,1,'Horse \rightarrow','FontSize',15); 
        % Add crossover info
            if joint == 2
                text(x0-.03,y0-.61,num2str(round(x0,2)),'FontSize',12)
            else
                text(x0-.07,y0-.61,num2str(round(x0,2)),'FontSize',12)
            end
        title([jointNames{joint},' - ',phaseString],'FontSize',16)
    end
%% Plot ALL TORQUES for a given length scale
figure('Position',[962,2,958,994]);
lengthNum = 9; range2plot = rangeMat(lengthNum,1):rangeMat(lengthNum,3); obj = objCell{lengthNum}; mid = rangeMat(lengthNum,2);
subplot(4,1,1)
    plot(obj.theta_motion_time(range2plot),jointmotion{lengthNum}(range2plot,:).*(180/pi)+[98.4373 102.226 116.2473],'LineWidth',3)
    title('Joint Motion','FontSize',16); ylabel('Joint Angle (deg)','FontSize',14); xlabel('Time (s)','FontSize',14); xlim([obj.theta_motion_time(range2plot(1)),obj.theta_motion_time(range2plot(end))])
    legend({'Hip';'Knee';'Ankle'}); xline(obj.theta_motion_time(mid),'HandleVisibility','off')
subplot(4,1,2)
    plot(obj.theta_motion_time(range2plot),torques_muscle{lengthNum}(range2plot,:),'LineWidth',3);
    title('Passive Muscle Torque','FontSize',16); ylabel('Torque (Nm)','FontSize',14); xlabel('Time (s)','FontSize',14); xlim([obj.theta_motion_time(range2plot(1)),obj.theta_motion_time(range2plot(end))])
    legend({'Hip';'Knee';'Ankle'});xline(obj.theta_motion_time(mid),'HandleVisibility','off')
subplot(4,1,3)
    plot(obj.theta_motion_time(range2plot),-torques_inertial{lengthNum}(range2plot,:),'LineWidth',3);
    title('Inertial Torque','FontSize',16); ylabel('Torque (Nm)','FontSize',14); xlabel('Time (s)','FontSize',14); xlim([obj.theta_motion_time(range2plot(1)),obj.theta_motion_time(range2plot(end))])
    legend({'Hip';'Knee';'Ankle'});xline(obj.theta_motion_time(mid),'HandleVisibility','off')
subplot(4,1,4)
    plot(obj.theta_motion_time(range2plot),-torques_active{lengthNum}(range2plot,:),'LineWidth',3);
    title('Total Joint Torque','FontSize',16); ylabel('Torque (Nm)','FontSize',14); xlabel('Time (s)','FontSize',14); xlim([obj.theta_motion_time(range2plot(1)),obj.theta_motion_time(range2plot(end))])
    legend({'Hip';'Knee';'Ankle'});xline(obj.theta_motion_time(mid),'HandleVisibility','off')
%% Plot Swing Phase Forces for a Zone
figure('Position',[1,1,1920,1003]); phaseInds = 654:1034; zoneNum = 3; lengthVals2plot = [1,5,8,11];
%legender = {'Rectus Femoris';'Vastus Lateralis';'Vastus Intermedius';'Vastus Medialis'};
nameTrim = @(x) x(4:end);
legender = cellfun(nameTrim,mDataZone{1}{zoneNum,2},'UniformOutput',false);
for ii = 1:length(lengthVals2plot)
    subplot(length(lengthVals2plot),1,ii)
    lengthNum = lengthVals2plot(ii);  
    timeVec = obj.theta_motion_time(phaseInds)-obj.theta_motion_time(phaseInds(1));
    swingpercVec = linspace(1,100,length(mDataZone{lengthNum}{zoneNum,1}(phaseInds,:)));
    dataVec = mDataZone{lengthNum}{zoneNum}(phaseInds,:)./max(mDataZone{lengthNum}{zoneNum}(phaseInds,:));
    plot(swingpercVec,dataVec,'LineWidth',3); grid on; set(gca,'FontSize',12)
    %legend(mDataZone{lengthNum}{zoneNum,2},'Interpreter','none','Location','eastoutside','FontSize',12); 
    legend(legender,'Interpreter','none','Location','eastoutside','FontSize',12); 
    title([num2str(round(lengthVals(lengthNum),2)),'x Rat Size'],'FontSize',18)
    ylabel({'Normalized';'Muscle Force'},'FontSize',16); xlabel('% Swing','FontSize',16); 
%     ylim([0 max(mDataZone{lengthNum}{zoneNum,1}(swingInds,:),[],'all')]); 
    ylim([0 1]);
    xlim([min(swingpercVec) max(swingpercVec)])
end
%% Plot Surface Figure for ONE Muscle
    dataLen = length(objCell{1}.theta_motion_time); stepInds = [654,1034,1676]; range2plot = 654:1676;
    mNum = 12;  amMat = zeros(numLens,length(range2plot)); yScaleLabels = 10;
    %figure
    for ii = 1:numLens
        %amMat(ii,:) = am_all{ii}(range2plot,mNum)./max(am_all{ii}(range2plot,mNum));
        if max(force_mn{ii}(range2plot,mNum)) == 0
            amMat(ii,:) = 0.*force_mn{ii}(range2plot,mNum);
        else
            amMat(ii,:) = force_mn{ii}(range2plot,mNum)./max(force_mn{ii}(range2plot,mNum));
        end
        %plot3(objCell{1}.theta_motion_time(range2plot),ii.*ones(1,length(am_all{ii}(range2plot,mNum))),am_all{ii}(range2plot,mNum)./max(am_all{ii}(range2plot,mNum))); hold on
    end
    
    % Make the activation matrix more "dense" with a beef factor which duplicates each activation waveform a desired number of times.
    bf = 3; temp = [];
    for ii = 1:numLens
        temp(bf*(ii-1)+1:bf*ii,:,1) = repmat(amMat(ii,:),[bf,1]);
    end
    amMat = temp(:,:,1);
    
    timeVec = objCell{1}.theta_motion_time(range2plot)-objCell{1}.theta_motion_time(range2plot(1));
    strideVec = linspace(0,100,length(timeVec));
    strideMat = repmat(strideVec,[numLens*bf, 1]);
    scaleMat = repmat(linspace(1,numLens*bf,numLens*bf)',[1 size(amMat,2)]);
    
    % For defining a custom color map between two colors
    color1 = [1,1,1]; color2 = [0,0,0];
    r = linspace(color1(1),color2(1),500)'; g = linspace(color1(2),color2(2),500)'; b = linspace(color1(3),color2(3),500)';  cm = [r,g,b];
    
    figure('Position',[1172,352,750,644]);
        surf(strideMat,scaleMat,amMat,'EdgeAlpha',0);view([0 90]); hold on; colormap parula %colormap(cm)
%         yTickVals = (bf:bf:numLens*bf)-floor(bf/2); yTickLabels = num2cell(linspace(lengthVals(1),lengthVals(end),length(yTickVals)));
        yTickLabelVals = round(lengthVals(floor(linspace(1,numLens,yScaleLabels))),2);
        yTickVals = floor(linspace(ceil(bf/2),numLens*bf-floor(bf/2),yScaleLabels)); yTickLabels = num2cell(yTickLabelVals);
        set(gca,'FontSize',14,'YTick',yTickVals,'YTickLabels',yTickLabels)
        % Add labels
            [~,ratInd] = min(abs(1-lengthVals)); [~,catInd] = min(abs(2-lengthVals)); [~,horseInd] = min(abs(12-lengthVals));
            text(101,bf*ratInd,1,'\leftarrow Rat','FontSize',15); yline(bf*ratInd,'k-.','LineWidth',1,'HandleVisibility','off','Alpha',.2)
            text(101,bf*catInd,1,'\leftarrow Cat','FontSize',15); yline(bf*catInd,'k-.','LineWidth',1,'HandleVisibility','off','Alpha',.2)
            text(101,bf*horseInd,1,'\leftarrow Horse','FontSize',15); yline(bf*horseInd,'k-.','LineWidth',1,'HandleVisibility','off','Alpha',.2)
            xline(strideVec(stepInds(2)-stepInds(1)),'c--','LineWidth',3,'HandleVisibility','off')
        xlim([0 max(strideVec)]); ylim([1 numLens*bf]); pbaspect([1 1 1]);
        ylabel('x Size Relative to Rat'); xlabel('Stride (%)') 
        title([objCell{1}.musc_obj{mNum}.muscle_name(4:end),' Activation v. Size Scale'],'FontSize',16)
%% Plot Surface Figure for a muscle GROUP
    zoneNum = 2;  yScaleLabels = 5; 
    stepInds = [beg,mid,ennd]; range2plot = beg:ennd; count  = 1; amMat = zeros(numLens,length(range2plot));  
    m2plot = find(cell2mat(muscZones(:,2))==zoneNum);
    
    zoneNames = {'HipAdductors';'HipFlexors';'KneeExtensors';'KneeFlexors';'AnklePlantarflexors';'AnkleDorsiflexors'};
    savePath = ['G:\My Drive\Rat\MeetingFiles\Meeting_20210712\SurfPlots\'];
    if exist(savePath, 'dir') ~= 7
       mkdir(savePath)
    end
    
    aa = figure('Position',[1,1,1920,1003],'Name',zoneNames{zoneNum});
    for kk = 1:length(m2plot)
        mNum = m2plot(kk);
        for ii = 1:numLens
            %amMat(ii,:) = am_all{ii}(range2plot,mNum)./max(am_all{ii}(range2plot,mNum));
            if max(force_mn{ii}(range2plot,mNum)) == 0
                amMat(ii,:) = 0.*force_mn{ii}(range2plot,mNum);
            else
                amMat(ii,:) = force_mn{ii}(range2plot,mNum)./max(force_mn{ii}(range2plot,mNum));
                %amMat(ii,:) = am_all{ii}(range2plot,mNum)./max(am_all{ii}(range2plot,mNum));
            end
            %plot3(objCell{1}.theta_motion_time(range2plot),ii.*ones(1,length(am_all{ii}(range2plot,mNum))),am_all{ii}(range2plot,mNum)./max(am_all{ii}(range2plot,mNum))); hold on
        end

        % Make the activation matrix more "dense" with a beef factor which duplicates each activation waveform a desired number of times.
        bf = 1; temp = [];
        for ii = 1:numLens
            temp(bf*(ii-1)+1:bf*ii,:,1) = repmat(amMat(ii,:),[bf,1]);
        end
        amMat = temp(:,:,1);

        timeVec = objCell{1}.theta_motion_time(range2plot)-objCell{1}.theta_motion_time(range2plot(1));
        strideVec = linspace(0,100,length(timeVec));
        strideMat = repmat(strideVec,[numLens*bf, 1]);
        scaleMat = repmat(linspace(1,numLens*bf,numLens*bf)',[1 size(amMat,2)]);

%         % For defining a custom color map between two colors
%         color1 = [1,1,1]; color2 = [0,0,0];
%         r = linspace(color1(1),color2(1),500)'; g = linspace(color1(2),color2(2),500)'; b = linspace(color1(3),color2(3),500)';  cm = [r,g,b];

        subplot(ceil(length(m2plot)/4),4,count);
            surf(strideMat,scaleMat,amMat,'EdgeAlpha',0);view([0 90]); hold on; colormap parula %colormap(cm)
            %set(gca,'yscale','log')
            %pPic = pcolor(strideMat,scaleMat,amMat); pPic.EdgeAlpha = 0;
            % Create accurate values for the Y Ticks
                yTickLabelVals = round(lengthVals(floor(linspace(1,numLens,yScaleLabels))),2);
                yTickVals = floor(linspace(ceil(bf/2),numLens*bf-floor(bf/2),yScaleLabels)); yTickLabels = num2cell(yTickLabelVals);
                set(gca,'FontSize',14,'YTick',yTickVals,'YTickLabels',yTickLabels)
            % Add labels
                [~,ratInd] = min(abs(1-lengthVals)); [~,catInd] = min(abs(2-lengthVals)); [~,horseInd] = min(abs(12-lengthVals)); [~,goatInd] = min(abs(6.46-lengthVals));
                xline(strideVec(stepInds(2)-stepInds(1)),'c--','LineWidth',3,'HandleVisibility','off')
                yline(bf*ratInd,'c-.','LineWidth',1,'HandleVisibility','off','Alpha',.5)
                yline(bf*catInd,'m-.','LineWidth',1,'HandleVisibility','off','Alpha',.5)
                yline(bf*horseInd,'y-.','LineWidth',1,'HandleVisibility','off','Alpha',.5)
                yline(bf*goatInd,'k-.','LineWidth',1,'HandleVisibility','off','Alpha',.5)
            if mod(kk,4) == 0 || (length(m2plot) < 4 && kk == length(m2plot))
                text(101,bf*ratInd+.2,1,'\leftarrow Rat','FontSize',15,'Color','c'); 
                text(101,bf*catInd+.2,1,'\leftarrow Cat','FontSize',15,'Color','m'); 
                text(101,bf*horseInd+.2,1,'\leftarrow Horse','FontSize',15,'Color',[145, 140, 39]./255);
                text(101,bf*goatInd+.2,1,'\leftarrow Goat','FontSize',15,'Color','k');
            end
            xlim([0 max(strideVec)]); ylim([1 numLens*bf]); pbaspect([1 1 1]);
            ylabel('x Size Relative to Rat'); xlabel('Stride (%)'); pbaspect([1 1 1])
            title([objCell{1}.musc_obj{mNum}.muscle_name(4:end)],'FontSize',16)
        count = count + 1;
    end
    saveas(aa,[savePath,zoneNames{zoneNum},'.png']) 
%% Plot STICK INSECT, RAT, HORSE example r2 figure
len = 11; joint = 1; rchLens = [8,13,22]; rchLabels = {'Stick Insect Sized';'Rat Sized';'Horse Sized'};
h = figure('Position',[7,427,1913,552]); 
timeVec = linspace(0,100,length(timeVec));
for ii = 1:3
    subplot(1,3,ii)
    beg = rangeMat(rchLens(ii),1); mid = rangeMat(rchLens(ii),2); ennd = rangeMat(rchLens(ii),3);
        timeVec = objCell{rchLens(ii)}.theta_motion_time(beg:ennd)-objCell{rchLens(ii)}.theta_motion_time(beg);
    timeVec = linspace(0,100,length(timeVec));
    ta = -torques_active{rchLens(ii)}(beg:ennd,joint); tv = -torques_muscle{rchLens(ii)}(beg:ennd,joint); ti = -torques_inertial{rchLens(ii)}(beg:ennd,joint);
    plot(timeVec,ta,'k','LineWidth',3); hold on;
    plot(timeVec,tv,'b','LineWidth',3);
    plot(timeVec,ti,'r','LineWidth',3);
    vCorr = corr(tv,ta);
    iCorr = corr(ti,ta);
    text(.7816,.1129,0,['r^2 = ',num2str(round(vCorr,2))],'FontSize',15,'Units','normalized','Color','blue'); 
    text(.7816,.0414,0,['r^2 = ',num2str(round(iCorr,2))],'FontSize',15,'Units','normalized','Color','red');
    xlim([0,max(timeVec)]); ylim([1.6*min([ta,ti,tv],[],'all'),1.1*max([ta,ti,tv],[],'all')])
    pbaspect([1,1,1]);
    ylabel('Torque (Nm)','FontSize',12)
    xlabel('Stride (%)','FontSize',12)
    title(rchLabels{ii},'FontSize',16)
    legend({'Active';'Viscoelastic';'Inertial'},'Location','southwest','FontSize',10)
end
saveas(h,['G:\My Drive\Rat\Swing Paper\r2_example'],'png')
%% SURF Torque Plots, 3 JOINTS
    aa = figure('Position',[1,1,1920,1003],'Name','Torques Across Scales'); range2plot = beg:ennd;yScaleLabels = 10; jointNames = {'Hip';'Knee';'Ankle'};
    for kk = 1:3
        tMat = torques_act{kk}';
        % Make the activation matrix more "dense" with a beef factor which duplicates each activation waveform a desired number of times.
        bf = 3; temp = [];
        for ii = 1:numLens
            temp(bf*(ii-1)+1:bf*ii,:,1) = repmat(tMat(ii,:),[bf,1]);
        end
        tMat = temp(:,:,1);

        %timeVec = objCell{1}.theta_motion_time(range2plot)-objCell{1}.theta_motion_time(range2plot(1));
        %strideVec = linspace(0,100,length(timeVec));
        strideVec = linspace(0,100,100);
        strideMat = repmat(strideVec,[numLens*bf, 1]);
        scaleMat = repmat(linspace(1,numLens*bf,numLens*bf)',[1 size(tMat,2)]);

        % For defining a custom color map between two colors
        color1 = [0,0,1]; color2 = [1,0,0];
        r = linspace(color1(1),color2(1),500)'; g = linspace(color1(2),color2(2),500)'; b = linspace(color1(3),color2(3),500)';  cm = [r,g,b];
        
        % For defining a custom color map between three colors
%         color1 = [0,0,1]; color2 = [1,1,1]; color3 = [1,0,0];
%         r = [linspace(color1(1),color2(1),250)';linspace(color2(1),color3(1),250)']; 
%         g = [linspace(color1(2),color2(2),250)';linspace(color2(2),color3(2),250)']; 
%         b = [linspace(color1(3),color2(3),250)';linspace(color2(3),color3(3),250)']; 
%         cm = [r,g,b];

        subplot(1,3,kk);
            surf(strideMat,scaleMat,tMat,'EdgeAlpha',0);view([0 90]); hold on; colormap(cm) %colormap parula %colormap(cm)
            %set(gca,'yscale','log')
            %pPic = pcolor(strideMat,scaleMat,amMat); pPic.EdgeAlpha = 0;
            % Create accurate values for the Y Ticks
                yTickLabelVals = round(lengthVals(floor(linspace(1,numLens,yScaleLabels))),2);
                yTickVals = floor(linspace(ceil(bf/2),numLens*bf-floor(bf/2),yScaleLabels)); 
                yTickLabels = num2cell(yTickLabelVals);
                set(gca,'FontSize',14,'YTick',yTickVals,'YTickLabels',yTickLabels)
            % Add labels
                [~,ratInd] = min(abs(1-lengthVals)); [~,catInd] = min(abs(2.66-lengthVals)); 
                [~,stickInd] = min(abs(.375-lengthVals)); [~,horseInd] = min(abs(12-lengthVals));
                xline(strideVec(37),'c--','LineWidth',3,'HandleVisibility','off')
                yline(bf*ratInd,'c-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*catInd,'m-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*horseInd,'y-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*stickInd,'w-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
            if mod(kk,3) == 0
                text(101,bf*ratInd,1,'\leftarrow Rat','FontSize',15); 
                text(101,bf*catInd,1,'\leftarrow Cat','FontSize',15); 
                text(101,bf*horseInd,1,'\leftarrow Horse','FontSize',15);
                text(101,bf*stickInd,1,'\leftarrow Stick Insect','FontSize',15);
            end
            colorbar('Location','southoutside','Ticks',[-1,1],'TickLabels',{'Flexion','Extension'})
            xlim([0 max(strideVec)]); ylim([1 numLens*bf]); pbaspect([1 1 1]);
            ylabel('x Size Relative to Rat'); xlabel('Stride (%)'); pbaspect([1 1 1])
            title([jointNames{kk},'- Active Joint Torque'],'FontSize',16)
    end
    %saveas(aa,[savePath,zoneNames{zoneNum},'.png'])
%% HIP SURF Torque Plot
    aa = figure('Position',[1,1,1920,1003],'Name','Torques Across Scales'); range2plot = beg:ennd;yScaleLabels = 10; jointNames = {'Hip';'Knee';'Ankle'};
    jointNum = 1;
    titleSize = 20; labelSize = 18; inWordSize = 16; lineWidth = 3; tickSize = 16;
    tMat = torques_act{jointNum}';
        % Make the activation matrix more "dense" with a beef factor which duplicates each activation waveform a desired number of times.
        bf = 3; temp = [];
        for ii = 1:numLens
            temp(bf*(ii-1)+1:bf*ii,:,1) = repmat(tMat(ii,:),[bf,1]);
        end
        tMat = temp(:,:,1);

        %timeVec = objCell{1}.theta_motion_time(range2plot)-objCell{1}.theta_motion_time(range2plot(1));
        %strideVec = linspace(0,100,length(timeVec));
        strideVec = linspace(0,100,100);
        strideMat = repmat(strideVec,[numLens*bf, 1]);
        scaleMat = repmat(linspace(1,numLens*bf,numLens*bf)',[1 size(tMat,2)]);

        % For defining a custom color map between two colors
        color1 = [0,0,1]; color2 = [1,0,0];
        r = linspace(color1(1),color2(1),500)'; g = linspace(color1(2),color2(2),500)'; b = linspace(color1(3),color2(3),500)';  cm = [r,g,b];
        cm = flipud(cool(500));
        
        % For defining a custom color map between three colors
%         color1 = [0,0,1]; color2 = [1,1,1]; color3 = [1,0,0];
%         r = [linspace(color1(1),color2(1),250)';linspace(color2(1),color3(1),250)']; 
%         g = [linspace(color1(2),color2(2),250)';linspace(color2(2),color3(2),250)']; 
%         b = [linspace(color1(3),color2(3),250)';linspace(color2(3),color3(3),250)']; 
%         cm = [r,g,b];

        subplot(1,2,1)
        torqueVec = -torques_active{13}(rangeMat(13,1):rangeMat(13,3),1); torqueVec = torqueVec./max(torqueVec);timeVec = linspace(0,100,length(torqueVec));
            plot(timeVec,torqueVec,'LineWidth',3)
            set(gca,'FontSize',14,'Position',[0.132,0.168,0.335,0.815])
            ylim([-1 1])
            xline(37,'r--','LineWidth',2.5); yline(0);
            set(gca,'YTick',[-1,0,1],'YTickLabels',{'FLX';'0';'EXT'},'FontSize',tickSize)
            set(gca,'XTick',0:20:100,'FontSize',tickSize)
            ylabel('Normalized Joint Torque','FontSize',labelSize); xlabel('Stride (%)','FontSize',labelSize) 
            title('Active Hip Torque - Rat','FontSize',titleSize)
            text(14.5,-0.92,'Swing','FontSize',16); text(64.5,-.92,'Stance','FontSize',inWordSize);
            text(-.1553,1.0482,1,'A','Units','normalized','FontSize',32)
            pbaspect([1,1,1])
        subplot(1,2,2);
            surf(strideMat,scaleMat,tMat,'EdgeAlpha',0);view([0 90]); hold on; colormap(cm) %colormap parula %colormap(cm)
            %set(gca,'yscale','log')
            %pPic = pcolor(strideMat,scaleMat,amMat); pPic.EdgeAlpha = 0;
            % Create accurate values for the Y Ticks
                yTickLabelVals = round(lengthVals(floor(linspace(1,numLens,yScaleLabels))),2);
                yTickVals = floor(linspace(ceil(bf/2),numLens*bf-floor(bf/2),yScaleLabels)); 
                yTickLabels = num2cell(yTickLabelVals);
                set(gca,'YTick',yTickVals,'YTickLabels',yTickLabels,'FontSize',tickSize)
            % Add labels
                [~,ratInd] = min(abs(1-lengthVals)); [~,catInd] = min(abs(2.66-lengthVals)); 
                [~,stickInd] = min(abs(.375-lengthVals)); [~,horseInd] = min(abs(12-lengthVals));
                xline(strideVec(37),'r--','LineWidth',3,'HandleVisibility','off')
                yline(bf*ratInd,'b-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*catInd,'g-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*horseInd,'y-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                yline(bf*stickInd,'w-.','LineWidth',2.5,'HandleVisibility','off','Alpha',.4)
                % Add text
                text(101,bf*ratInd,1,'\leftarrow Rat','FontSize',tickSize); 
                text(101,bf*catInd,1,'\leftarrow Cat','FontSize',tickSize); 
                text(101,bf*horseInd,1,'\leftarrow Horse','FontSize',tickSize);
                text(101,bf*stickInd,1,'\leftarrow Stick Insect','FontSize',tickSize);
            colorbar('Location','southoutside','Ticks',[-1,1],'TickLabels',{'Flexion','Extension'},'FontSize',labelSize)
            xlim([0 max(strideVec)]); ylim([1 numLens*bf]); pbaspect([1 1 1]);
            ylabel('x Size Relative to Rat','FontSize',labelSize); xlabel('Stride (%)','FontSize',labelSize);
            text(-.1258,1.048,1,'B','Units','normalized','FontSize',32)
            title(['Active Hip Torque - All Scales'],'FontSize',titleSize)
%% Calculate Swing Data
stepInds = [654,1034,1676]; range2plot = stepInds(1):stepInds(3); ft = []; legender = cell(1,1); showStance = 0; clear sw

    if showStance
        phaseInds = [floor((1/3)*(stepInds(3)-stepInds(2))),...
                     floor((2/3)*(stepInds(3)-stepInds(2))),...
                     (stepInds(3)-stepInds(2))];
        strideName = 'Stance';
    else
        phaseInds = [floor((1/3)*(stepInds(2)-stepInds(1))),...
                     floor((2/3)*(stepInds(2)-stepInds(1))),...
                     (stepInds(2)-stepInds(1))];
        strideName = 'Swing';
    end

for joint = 1
    if showStance
        ft = torques_act{joint}(stepInds(2):stepInds(3),:);
    else
        ft = torques_act{joint}(stepInds(1):stepInds(2),:);
    end
    sw(:,1) = mean(ft(1:phaseInds(1),:));
    sw(:,2) = mean(ft(phaseInds(1):phaseInds(2),:));
    sw(:,3) = mean(ft(phaseInds(2):phaseInds(3),:));
end

% Define colormap
cm = hsv(length(lengthVals));

lens2plot = [5,11,20,25];
lens2plot = 11:2:25;
lens2plot = 1:numLens;
figure('Position',[2,2,958,994],'Name','Torques Across Scales'); count = 1;
bPlot = bar(sw(lens2plot,:)','FaceAlpha',1,'EdgeAlpha',0,'FaceColor','flat');hold on;
    % Set scale colors
    for ii = 1:length(lens2plot)
       bPlot(ii).CData = cm(lens2plot(ii),:); 
    end
    % Tweak tick info on axes
    set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
    set(gca,'YTick',-1:.5:1,'YTickLabel',{'Flexion  -1';'-.5';'0';'.5';'Extension  1'})
    % Add labels and title
    legend(cellfun(@num2str,num2cell(round(lengthVals(lens2plot),2)),'UniformOutput',false),'Location','eastoutside'); ylim([-1, 1]); grid on
    ylabel('Normalized Joint Torque','FontSize',18); xlabel('Stride Phase','FontSize',18)
    title(['Active Joint Torques of the ',objCell{1}.joint_obj{joint}.name(4:end)],'FontSize',20)
%% Normalize a torque matrix
function outMat = normalizeTorque(torque_temp,trans,allTogether)
    if allTogether
        torque_temp = interp1(1:length(torque_temp),torque_temp,linspace(1,length(torque_temp),100));
        outMat = torque_temp./max(abs(torque_temp));
        return
    end
    if size(torque_temp,2) == 1
        oneHundredizer = @(inMat,trans) [interp1(1:length(inMat(1:trans,:)),inMat(1:trans,:),linspace(1,length(inMat(1:trans,:)),37)),...
                                         interp1(1:length(inMat(trans:end,:)),inMat(trans:end,:),linspace(1,length(inMat(trans:end,:)),63))]';
    else
        oneHundredizer = @(inMat,trans) [interp1(1:length(inMat(1:trans,:)),inMat(1:trans,:),linspace(1,length(inMat(1:trans,:)),37));...
                                         interp1(1:length(inMat(trans:end,:)),inMat(trans:end,:),linspace(1,length(inMat(trans:end,:)),63))];
    end
    torque_temp = oneHundredizer(torque_temp,trans);  
    if 0
        swingNorm = -1 + 2.*(torque_temp(1:37,:) - min(torque_temp(1:37,:)))./(max(torque_temp(1:37,:)) - min(torque_temp(1:37,:)));
        stanceNorm = -1 + 2.*(torque_temp(38:end,:) - min(torque_temp(38:end,:)))./(max(torque_temp(38:end,:)) - min(torque_temp(38:end,:)));
        outMat = [swingNorm;stanceNorm];
    else
        swingNorm = torque_temp(1:37,:)./max(abs(torque_temp(1:37,:)));
        stanceNorm = torque_temp(38:end,:)./max(abs(torque_temp(38:end,:)));
        outMat = [swingNorm;stanceNorm];
        %outMat = torque_temp./max(abs(torque_temp));
    end
end
%% Find_Crossover_Point
function [outX,outY] = find_crossover_point(corr1,corr2,lengthVals)
    % Scan across the two correlation vectors and find the point at which they intersect
    aa = corr1; bb = corr2;    
    for ii = 2:length(aa)
        if aa(ii) < bb(ii) && aa(ii-1) > bb(ii-1) % if we just passed the crossover point, find the intersection
            m2 = (bb(ii)-bb(ii-1))/(lengthVals(ii)-lengthVals(ii-1));
            x = lengthVals(ii-1);
            y2 = bb(ii-1);
            m1 = (aa(ii)-aa(ii-1))/(lengthVals(ii)-lengthVals(ii-1));
            y1 = aa(ii-1);
            outX = ((m2-m1)*x+y1-y2)/(m2-m1);
            outY = (m2*y1-m1*y2)/(m2-m1);
        end
    end
end
%% Range_finder: Output Indices for Range of Motion
function range2plot = range_finder(range,len,rangeMat)
    beg = rangeMat(len,1); mid = rangeMat(len,2); ennd = rangeMat(len,3);
    switch range
        case 1
            endInd = beg + floor((mid-beg)/3);
            range2plot = beg:endInd;
        case 2
            endInd1 = beg + floor((mid-beg)/3);
            endInd2 = beg + floor(2*(mid-beg)/3);
            range2plot = endInd1:endInd2;
        case 3
            endInd = beg + floor(2*(mid-beg)/3);
            range2plot = endInd:mid;
        case 4
            range2plot = rangeMat(len,1):rangeMat(len,2);
        case 5
            endInd = mid + floor((ennd-mid)/3);
            range2plot = mid:endInd;
        case 6
        case 7

        case 8
            range2plot = rangeMat(len,2):rangeMat(len,3);
    end
end