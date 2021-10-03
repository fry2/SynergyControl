curVars = whos;
if ~any(contains({curVars.name},'AllDat'),'all')
    load('C:\Users\fry16\OneDrive\Documents\NW_data_local\AllDat.mat')
end

if ~any(contains({curVars.name},'outScale'),'all')
    load([pwd,'\scalingOpt.mat'],'outScale')
end

extScaling = outScale;

trialNum = 2;

trialData = AllDat{trialNum};
emgData = trialData.emg_env;
jointData = trialData.joints;
emgNames = trialData.emg_names;

if size(emgData,3) == size(jointData,3)
    numSteps = size(jointData,3);
else
    error('Number of steps don''t match between provided EMG signals and the joint motion')
end

timeVec = 0:.54e-3:10.01;
simPath = 'G:\My Drive\Rat\SynergyControl\Data\EMG_walking_baseSim.asim';
%outScale = zeros(10,12);

tstart = tic;
for ii = 2%numSteps
    if sum(extScaling(ii,:)) ~= 0
        startScale = extScaling(ii,:);
    else
        startScale = [];
    end
    emgStep = emgData(:,:,ii);
    if ~(sum(isnan(mean(emgStep)))>3)
        tstartStep = tic;
        if any(isnan(mean(emgStep)))
            emgStep(isnan(emgStep)) = 0;
        end
        inEMGdata{1} = emgStep;
        inEMGdata{2} = emgNames;
        inEMGdata = order_emg_for_sim(simPath,inEMGdata);
        jointMotion = jointData(:,[1,3,4],ii);
        %jointMotion = [jointMotion;jointMotion;jointMotion];
        jointBig = interp1(1:length(jointMotion),jointMotion,linspace(1,length(jointMotion),length(timeVec)));
        [outScale(ii,:),scaleVals(ii)] = emg_walking_optimization(inEMGdata,jointBig,startScale);
        telapsedStep = toc(tstartStep);
        disp(['Step ',num2str(ii),' ',num2str(telapsedStep),' seconds.'])
    end
end
telapsed = toc(tstart);
disp([num2str(telapsed/60),' minutes for full optimization.'])

save([pwd,'\EMG_control\scalingOpt.mat'],'outScale')

%% 
inSimPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking_Standalone.asim';
inScale = outScale(2,:);
jointMotion = emg_to_walking_simproj(inSimPath,inEMGdata,inScale);

%% some plotting, can delete
figure;subplot(3,1,1); plot(linspace(0,10,length(inEMGdata{1})),inEMGdata{1},'LineWidth',1.5);title('Activation');xlabel('Time (s)')
subplot(3,1,2);plot(timeVec,jointBig,'LineWidth',2);xlim([0 10]);title('Desired Joint Motion');xlabel('Time (s)')
subplot(3,1,3);plot(timeVec(1:length(jointMotion)),jointMotion.*(180/pi)+[98.4373 102.226 116.2473],'LineWidth',2);xlim([0 10]);title('Resulting Joint Motion');xlabel('Time (s)')