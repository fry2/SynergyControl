inSim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_neutral.asim";
% an insim with constant joint positions held at a neutral position: [0.089757, 0.25314, 0.38158]
sdata = processSimData(inSim);
lenInd = find(contains({sdata.name},'KeyMuscleLen'));
lendata = sdata(lenInd).data(1000:end-20,:);
muscnames = sdata(lenInd).data_cols;
Lwidth_func = @(perc1,Lr) Lr.*sqrt(.25./(1-perc1)); % Based on percentage value when L = 1.5Lr (ex. perc = .7)


%%
% preallocate
numIts = 50; clear pTorques
lrs = zeros(38,numIts); lmins = lrs; lmaxs = lrs; lwidths = lrs; llim = .7; ulim = 1;

for ii = 1:numIts
    rng(ii)
    r = llim + (ulim-llim).*rand(1,38);
    lrs(:,ii) = (r.*mean(lendata))';
    lmins(:,ii) = .5.*lrs(:,ii);
    lmaxs(:,ii) = 1.5.*lrs(:,ii);
    lwidths(:,ii) = Lwidth_func(.7,lrs(:,ii));
end

    lrs(:,end+1) = mean(lendata)';
    lmins(:,end+1) = .5.*lrs(:,end);
    lmaxs(:,end+1) = 1.5.*lrs(:,end);
    lwidths(:,end+1) = Lwidth_func(.7,lrs(:,end));

%%

% Consider loading data: load([pwd,'/Data/lr_randomizer_sensitivty_data_10trials_p5_to_1p5.mat'],'pTorques','jm','locs','lrs','lmins','lmaxs','lwidths')
inText = importdata(inSim);
cStimNames = {'ConstantHip';'ConstantKnee';'ConstantxAnkle'}; count = 0;
for ii = 1:numIts+1
    tstart = tic;
    for muscNum = 1:38
        muscInd = find(contains(inText,['<Name>LH_',muscnames{muscNum},'</Name>'],'IgnoreCase',true));
        lrInd = find(contains(inText(muscInd:end),'RestingLength'),1,'first')+muscInd-1;
        lwInd = find(contains(inText(muscInd:end),'Lwidth'),1,'first')+muscInd-1;
        llInd = find(contains(inText(muscInd:end),'LowerLimit'),1,'first')+muscInd-1;
        ulInd = find(contains(inText(muscInd:end),'UpperLimit'),1,'first')+muscInd-1;
        inText{lrInd} = replaceBetween(inText{lrInd},'>','</',num2str(lrs(muscNum,ii)));
        inText{lwInd} = replaceBetween(inText{lwInd},'>','</',num2str(lwidths(muscNum,ii)));
        inText{llInd} = replaceBetween(inText{llInd},'>','</',num2str(lmins(muscNum,ii)));
        inText{ulInd} = replaceBetween(inText{ulInd},'>','</',num2str(lmaxs(muscNum,ii)));
    end
    for stimNum = 1:3
        stimInd = find(contains(inText,cStimNames{stimNum}));
        inText{stimInd+2} = '<Enabled>False</Enabled>';
    end
    [~,jobID] = fileparts(tempname);
    tempPath = ['C:\Users\fry16\OneDrive\Documents\MATLAB\tempFolder\',jobID,'.asim'];
    fileID = fopen(tempPath,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);
    [~,~,~,~,~,objH] = muscle_to_joint_VE(tempPath,1);
    tempH = compute_passive_joint_torque(objH);
    [~,~,~,~,~,objK] = muscle_to_joint_VE(tempPath,2);
    tempK = compute_passive_joint_torque(objK);
    [~,~,~,~,~,objA] = muscle_to_joint_VE(tempPath,3);
    tempA = compute_passive_joint_torque(objA);
    pTorques(:,ii,1) = tempH(1,:)';
    pTorques(:,ii,2) = tempK(2,:)';
    pTorques(:,ii,3) = tempA(3,:)';
    %pTorqueDoc(:,:,ii) = [pTorqueH(1,:)',pTorqueK(2,:)',pTorqueA(3,:)'];
    count  = count + 1; telapsed = toc(tstart);
    disp(['Finished ',num2str(count),' out of ',num2str(numIts),'. (',num2str(telapsed),' s).'])
end
[~,locsH] = findpeaks(-objH.theta_motion(:,1));
[~,locsK] = findpeaks(objK.theta_motion(:,2));
[~,locsA] = findpeaks(-objA.theta_motion(:,3));
locs = [locsH,locsK,locsA]; jm = [objH.theta_motion(:,1),objK.theta_motion(:,2),objA.theta_motion(:,3)];
%% Plot the results, mean and std
figure('Position',[962,2,958,994]);

titles = {{['Passive Torque for ',num2str(size(pTorques,2)-1),' Random Resting Lengths, (',num2str(llim),'-',num2str(ulim),'):'];'Hip'};'Knee';'Ankle'};
%bnds = [-.2, .31;-.1 .05;-.02 .04];
for jNum = 1:3
    subplot(3,1,jNum)
    tProf = pTorques(locs(1,jNum):locs(2,jNum),1:end-1,jNum); tLen = length(locs(1,jNum):locs(2,jNum));
    stUp = (mean(tProf,2)+std(tProf,[],2))'; stDn = (mean(tProf,2)-std(tProf,[],2))';
    xVec = linspace(0,100,length(tProf)); 
    patch([xVec fliplr(xVec)], [stDn,fliplr(stUp)], 'k','FaceAlpha',.2,'EdgeAlpha',0); hold on
    plot(xVec,mean(tProf(:,1:end-1),2),'k','LineWidth',3);
    plot(xVec,pTorques(locs(1,jNum):locs(2,jNum),end,jNum),'r','LineWidth',2)
    plot(xVec,pTorques(locs(1,jNum):locs(2,jNum),1,jNum),'b','LineWidth',2)
    xlim([0 max(xVec)]); grid on;ylim([min(tProf,[],'all') max(tProf,[],'all')])
    title(titles{jNum},'FontSize',14); xlabel('Percent Stride (%)'); ylabel('Passive Torque (N-m)'); 
    legend({'Standard Deviation';'Mean';'Without Randomness';'Current Config'},'Location','southeast','FontSize',8)
end
%% Plot the results, all waves
figure('Position',[962,2,958,994]);

titles = {{['Passive Torque for ',num2str(size(pTorques,2)-1),' Random Resting Lengths, (',num2str(llim),'-',num2str(ulim),'):'];'Hip'};'Knee';'Ankle'};
%bnds = [-.2, .31;-.1 .05;-.02 .04];
for jNum = 1:3
    subplot(3,1,jNum)
    tProf = pTorques(locs(1,jNum):locs(2,jNum),1:end-1,jNum); tLen = length(locs(1,jNum):locs(2,jNum));
    stUp = (mean(tProf,2)+std(tProf,[],2))'; stDn = (mean(tProf,2)-std(tProf,[],2))';
    xVec = linspace(0,100,length(tProf));
    plot(xVec,tProf(:,1:end-1)); hold on
    plot(xVec,pTorques(locs(1,jNum):locs(2,jNum),end,jNum),'r','LineWidth',3)
    plot(xVec,pTorques(locs(1,jNum):locs(2,jNum),1,jNum),'b','LineWidth',2)
    xlim([0 max(xVec)]); grid on;ylim([min(tProf,[],'all') max(tProf,[],'all')])
    title(titles{jNum},'FontSize',14); xlabel('Percent Stride (%)'); ylabel('Passive Torque (N-m)'); 
end
%% which group of lrs is the closest to the mean line?
valDoc = [(1:50)'];
for jNum = 1:3
    tProf = pTorques(locs(1,jNum):locs(2,jNum),1:end-1,jNum); tLen = length(locs(1,jNum):locs(2,jNum));
    meanT = mean(tProf,2);
    %valDoc = [valDoc,sortrows([(1:50)',sum((tProf-meanT).^2)'],2)];
    valDoc(:,jNum+1) = sum((tProf-meanT).^2)';
end
valDoc(:,5) = sum(valDoc(:,2:4),2);
[minVal,minLoc] = min(valDoc(:,5));
disp(['Minimum value from lrs column ',num2str(valDoc(minLoc,1))])
for jNum = 1:3
    tProf = pTorques(locs(1,jNum):locs(2,jNum),1:end-1,jNum); meanT = mean(tProf,2);
    subplot(3,1,jNum)
    plot(meanT,'k'); hold on
    plot(tProf(:,valDoc(minLoc,1)),'b')
    plot(pTorques(locs(1,jNum):locs(2,jNum),1,jNum),'r')
    legend({'Mean';'New';'Current'})
end
%% Save these values into neutral_lengths
neutral_lengths = struct();
neutral_lengths.data_cols = {'Muscle Names','Lr','Lw','Lmin','Lmax','Perc'};
neutral_lengths.data = [muscnames',num2cell(lrs(:,minLoc)),num2cell(lwidths(:,minLoc)),num2cell(lmins(:,minLoc)),num2cell(lmaxs(:,minLoc)),num2cell(.7*ones(38,1))];
save([pwd,'\Data\neutral_lengths.mat'],'neutral_lengths')
disp('Data saved.')