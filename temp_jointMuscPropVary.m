tstart = tic;
simPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim";
%trialPath = 'G:\My Drive\Rat\SynergyControl\JointStiffnessAnalysis\HipStiffness\';
trialPath = 'C:\Users\fry16\OneDrive\Desktop\HipStiffness\';
delete([trialPath,'*.*'])

simContents = importdata(simPath);
%Hip
muscTags = {'tfl';'ab'};
muscNames = {'LH_TensorFasciaLatae';'LH_AdductorBrevis'};
musc1stimID = 'stTC1-1TFL';
musc2stimID = 'stTC1-7AB';
m1Steep_0 = 1382.67;
m2Steep_0 = 1681.89;
jointFricID = '<ID>028bd626-f964-4165-8e5b-e5f69937ba5e</ID>';
jointNum = 3; %as decided by the JointMotion.txt columns

%Knee
% muscTags = {'bfp';'rf'};
% musc1stimID = 'stTC1-2BFP';
% musc2stimID = 'stTC1-3RF';
% m1Steep_0 = 1273.72;
% m2Steep_0 = 1607.74;
% jointFricID = '<ID>00b42b7b-a001-4626-a69c-87c0beca56b3</ID>';
% jointNum = 2;

%Ankle
% muscTags = {'ta';'mg'};
% musc1stimID = 'stTC1-5TA';
% musc2stimID = 'stTC1-4MG';
% m1Steep_0 = 1251.94;
% m2Steep_0 = 2713.95;
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
stiffCoeffs = linspace(0,.75,numStiffs);
%stiffCoeffs = 0;

fun = @(inVec) objFunc(inVec,simContents,m1StimVals,m2StimVals,m1Ind,m2Ind,fricInd,trialPath,muscTags,jointNum,muscNames);
%     options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-10,'FiniteDifferenceStepSize',.01,'TolCon',1e-8,'Display','iter-detailed');
     options = optimoptions('fmincon','Algorithm','sqp','Display','iter-detailed',...
         'PlotFcn',{'optimplotx';'optimplotfval';'optimplotstepsize'});
    %options = optimoptions('fmincon','Algorithm','sqp','FiniteDifferenceType','central','UseParallel',true,'Display','iter-detailed');
    findIO = 0;
    if findIO
        numTests =  5;
        counter = 1;
        bestVal = 1000;
        while counter <= numTests
            testI = [m1Steep_0*rand(1,1) m2Steep_0*rand(1,1) rand(1,1)];
            testOut = fun(testI);
            if testOut < bestVal
                bestVal = testOut;
                i0 = testI;
            end
            counter = counter+1;
            if ~mod(counter,20)
                disp(num2str(counter))
            end
        end
        delete([trialPath,'*.*'])
    else
        i0 = [m1Steep_0,.5*m2Steep_0,1];
    end
    lb = zeros(1,3);
    ub = [1e5,1e5,1];
    outI = fmincon(fun,i0,[],[],[],[],lb,ub,[],options);
    [outVal,stiffCell] = fun(outI);

[m1plotStim,m2plotStim] = meshgrid(m1StimVals,m2StimVals);
zlims = [min(stiffCell,[],'all') max(stiffCell,[],'all')];

figure
surf(m1plotStim,m2plotStim,stiffCell');
pbaspect([1,1,1])
xlabel([muscTags{1},' Stim'])
ylabel([muscTags{2},' Stim'])
zlabel(['Stiffness ',num2str(50000)])
zlim(zlims);
view([45 40])

    telapsed = toc(tstart);
    disp(['Calculation time:',' ',num2str(telapsed),'s for ',num2str(length(resCell)),' calculations. (',num2str(telapsed/length(resCell)),'s per).'])
            
function [outVal,stiffCell] = objFunc(inVec,simContents,m1StimVals,m2StimVals,m1Ind,m2Ind,fricInd,trialPath,muscTags,jointNum,muscNames)
    
    steep1 = inVec(1);
    steep2 = inVec(2);
    stiffness = inVec(3);
    % muscName = '<Name>LH_TensorFasciaLatae</Name>';
    for tt = 1:length(muscNames)
        muscName = muscNames{tt};
        temp = find(contains(simContents,['<Name>',muscName,'</Name>']));
        muscSteepInd = find(contains(simContents(temp:end),'<C>'),1,'first')+temp-1;
        switch tt
            case 1
                newSteepVal = steep1;
            case 2
                newSteepVal = steep2;
        end
        simContents{muscSteepInd} = replaceBetween(simContents{muscSteepInd},'>','<',num2str(newSteepVal));
    end
    numStims = length(m1StimVals);

    simContents{fricInd} = replaceBetween(simContents{fricInd},'>','<',num2str(stiffness));
    for ii = 1:numStims
        m1Stim = m1StimVals(ii);
        for jj = 1:numStims
            m2Stim = m2StimVals(jj);
            simContents{m1Ind} = replaceBetween(simContents{m1Ind},'>','<',[num2str(m1Stim/10),'e-008']);
            simContents{m2Ind} = replaceBetween(simContents{m2Ind},'>','<',[num2str(m2Stim/10),'e-008']);
            txtInd = find(contains(simContents,'OutputFilename>Joint'));
            simContents{txtInd} = replaceBetween(simContents{txtInd},'ion','.txt',['_',num2str(m1Stim),'_',num2str(m2Stim),'_',num2str(stiffness)]);
            trialFileName = [trialPath,muscTags{1},num2str(m1Stim),'_',muscTags{2},num2str(m2Stim),'_stiff',num2str(stiffness),'.asim'];
            fileID = fopen(trialFileName,'w');
            fprintf(fileID,'%s\n',simContents{:});
            fclose(fileID);
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
    
    stiffCell = zeros(length(m1StimVals),length(m1StimVals));
    for ii = 1:length(resCell)
        fileString = string(resCell{ii}{1});
        m1Val = str2double((extractBetween(fileString,muscTags{1},'_')));
            [~,m1Ind] = min(abs(m1StimVals-m1Val));
        m2Val = str2double((extractBetween(fileString,muscTags{2},'_')));
            [~,m2Ind] = min(abs(m2StimVals-m2Val));
        stiffVal = str2double((extractBetween(fileString,'stiff','.asim')));
        stiffCell(m1Ind,m2Ind) = resCell{ii}{3};
        %stiffRes(m1Ind,m2Ind) = resCell{ii}{2};
    end
    
    jointLims = [-70 11;...
                 -77 3;...
                 -94 5]*(pi/180);
             
    jLim = jointLims(jointNum,:);
    
    resBnds = [min(stiffCell,[],'all') max(stiffCell,[],'all')];
    
    outVal = sum((resBnds-jLim).^2) + 20*max(gradient(stiffCell),[],'all') + 6*sum(gradient(stiffCell)./max(gradient(stiffCell),[],'all')<=.005,'all')/numel(stiffCell);
end
