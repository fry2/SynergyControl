tstart = tic;
simPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim";
%trialPath = 'G:\My Drive\Rat\SynergyControl\JointStiffnessAnalysis\HipStiffness\';
trialPath = 'C:\Users\fry16\OneDrive\Desktop\HipStiffness\';
delete([trialPath,'*.*'])

simContents = importdata(simPath);
%Hip
muscTags = {'tfl';'ab'};
musc1stimID = 'stTC1-1TFL';
musc2stimID = 'stTC1-7AB';
jointFricID = '<ID>028bd626-f964-4165-8e5b-e5f69937ba5e</ID>';
jointNum = 3; %as decided by the JointMotion.txt columns

%Knee
% muscTags = {'bfp';'rf'};
% musc1stimID = 'stTC1-2BFP';
% musc2stimID = 'stTC1-3RF';
% jointFricID = '<ID>00b42b7b-a001-4626-a69c-87c0beca56b3</ID>';
% jointNum = 2;

%Ankle
% muscTags = {'ta';'mg'};
% musc1stimID = 'stTC1-5TA';
% musc2stimID = 'stTC1-4MG';
% jointFricID = '<ID>064c80b2-4f34-4dac-9298-096da5037c32</ID>';
% jointNum = 1;

m1Ind = find(contains(simContents,musc1stimID))+13;
m2Ind = find(contains(simContents,musc2stimID))+13;
fricInd = find(contains(simContents,jointFricID))+4;
%For modifying the static friction ratio
%fricInd = find(contains(simContents,'<ID>028bd626-f964-4165-8e5b-e5f69937ba5e</ID>'))+7;

numStims = 5;
numStiffs = 5;

m1StimVals = linspace(0,20,numStims);
m2StimVals = linspace(0,20,numStims);
stiffCoeffs = linspace(0,1,numStiffs);
stiffCoeffs = .25;

for kk = 1:length(stiffCoeffs)
    fric = stiffCoeffs(kk);
    simContents{fricInd} = replaceBetween(simContents{fricInd},'>','<',num2str(fric));
    for ii = 1:numStims
        m1Stim = m1StimVals(ii);
        for jj = 1:numStims
            m2Stim = m2StimVals(jj);
            simContents{m1Ind} = replaceBetween(simContents{m1Ind},'>','<',[num2str(m1Stim/10),'e-008']);
            simContents{m2Ind} = replaceBetween(simContents{m2Ind},'>','<',[num2str(m2Stim/10),'e-008']);
            txtInd = find(contains(simContents,'OutputFilename>Joint'));
            simContents{txtInd} = replaceBetween(simContents{txtInd},'ion','.txt',['_',num2str(m1Stim),'_',num2str(m2Stim),'_',num2str(fric)]);
            trialFileName = [trialPath,muscTags{1},num2str(m1Stim),'_',muscTags{2},num2str(m2Stim),'_stiff',num2str(fric),'.asim'];
            fileID = fopen(trialFileName,'w');
            fprintf(fileID,'%s\n',simContents{:});
            fclose(fileID);
        end
    end
end
trialFiles = dir(trialPath);
trialFiles = {trialFiles(contains({trialFiles.name},'.asim')).name};
resCell = cell(size(trialFiles));
    %sour_folder = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin';
    sour_folder = 'C:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin';

% % Now do the actual simulations
    parfor i = 1:length(trialFiles)
        fileName = trialFiles{i};
        executable = ['"',sour_folder,'\AnimatSimulator" "',[trialPath,fileName],'"'];
%         [status, message] = system(executable);
        jsystem(executable,'noshell');
        m1Val = str2double((extractBetween(fileName,muscTags{1},'_')));
        m2Val = str2double((extractBetween(fileName,muscTags{2},'_')));
        stiffVal = str2double((extractBetween(fileName,'stiff','.asim')));
        ds = importdata([trialPath,['JointMotion_',num2str(m1Val),'_',num2str(m2Val),'_',num2str(stiffVal),'.txt']]);
        jointProfile = ds.data(:,3:5);
        resCell{i} = [{fileName},{jointProfile},{mean(jointProfile(2000:6000,jointNum))}];
    end

stiffCell = cell(1,length(stiffCoeffs));
stiffRes = stiffCell;
for ii = 1:length(resCell)
    fileString = string(resCell{ii}{1});
    m1Val = str2double((extractBetween(fileString,muscTags{1},'_')));
        [~,m1Ind] = min(abs(m1StimVals-m1Val));
    m2Val = str2double((extractBetween(fileString,muscTags{2},'_')));
        [~,m2Ind] = min(abs(m2StimVals-m2Val));
    stiffVal = str2double((extractBetween(fileString,'stiff','.asim')));
    [~,stiffBin] = min(abs(stiffCoeffs-stiffVal));
    temp = size(stiffCell{stiffBin},1);
    stiffCell{stiffBin}(m1Ind,m2Ind) = resCell{ii}{3};
    stiffRes{stiffBin}{m1Ind,m2Ind} = resCell{ii}{2};
end

[m1plotStim,m2plotStim] = meshgrid(m1StimVals,m2StimVals);
zlims = [min(cell2mat(stiffCell),[],'all') max(cell2mat(stiffCell),[],'all')];

for ii = 1:length(stiffCoeffs)
    subplot(1,length(stiffCoeffs),ii)
surf(m1plotStim,m2plotStim,stiffCell{ii}');
pbaspect([1,1,1])
xlabel([muscTags{1},' Stim'])
ylabel([muscTags{2},' Stim'])
zlabel(['Stiffness ',num2str(stiffCoeffs(ii))])
zlim(zlims);
view([45 40])
end

            telapsed = toc(tstart);
            disp(['Calculation time:',' ',num2str(telapsed),'s for ',num2str(length(resCell)),' calculations. (',num2str(telapsed/length(resCell)),'s per).'])

if length(stiffCoeffs)>1            
    for ii = 1:length(resCell)
        initPos(ii) = mean(resCell{ii}{2}(1:20,2));
    end
    for ii=1:length(stiffCell)
        coStimPos(ii) = stiffCell{ii}(end,end);
        stiffRange(ii) = max(stiffCell{ii},[],'all')-min(stiffCell{ii},[],'all');
    end
    [~,temp1] = min(abs(coStimPos-mean(initPos)));
    [~,temp2] = max(stiffRange);
    disp(['Stiffness where mutual full stimulation is closest to initial position: ',num2str(stiffCoeffs(temp1))])
    disp(['Stiffness with largest range: ',num2str(stiffCoeffs(temp2))])
%     figure; plot(stiffCoeffs,coStimPos)
%     hold on
%     plot(stiffCoeffs,mean(initPos).*ones(1,length(stiffCoeffs)))
end
            
% For interploating surfaces
% [tflplotStimBig,abplotStimBig] = meshgrid(1:.1:20,1:.1:20);
% stiffCellBig = griddata(tflplotStim,abplotStim,stiffCell{5},tflplotStimBig,abplotStimBig);
