neutralSim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_neutral.asim";
nobj = design_synergy(neutralSim);
bodyTorque = compute_body_torques(nobj,0);
%%
for ii = 1:3
    temp = compute_joint_moment_arms(nobj,ii,1);
    moment_arms(:,:,ii) = temp(1:38,:)'./1000;
end
%%
for ii = 1:38
    musc = nobj.musc_obj{ii}; ml = musc.muscle_length_profile; lr = musc.RestingLength; ks = musc.Kse; kp = musc.Kpe;
    k0(ii,1) = (ks*kp)/(ks+kp);
    for jj = 1:3
        Atot = max(ml-lr,0).*moment_arms(:,ii,jj);
        Aeq(jj,ii) = mean(Atot(3e3:15e3));
    end
    muscDat{ii,1} = musc.muscle_name;
    muscDat{ii,2} = k0(ii,1);
end
%%
beq = -mean(bodyTorque(3e3:15e3,:))';
%lb = max(1.*ones(38,1),k0);
lb = 10.*ones(38,1);

kfunc = @(x) sum((x-k0).^2);
nonlcon_wrap = @(x) k_compute_nonlcon_test(x,inText,muscDat);

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','sqp');

[x,fval] = fmincon(kfunc,k0,[],[],[],[],[],[],nonlcon_wrap,options);
muscDat(:,3) = num2cell(x);

%%
figure;bar(k0,'FaceAlpha',.5,'EdgeAlpha',0);hold on;bar(x,'FaceAlpha',.5,'EdgeAlpha',0);legend({'Old Vals';'New Vals'})

%% With new values, put it into a sim file and run it, does it work?
    inText = importdata(neutralSim);

    % Update the K values
    for ii = 1:length(muscDat)
        muscInd = find(contains(inText,['<Name>',muscDat{ii,1},'</Name>'],'IgnoreCase',true));
        ksInd = find(contains(inText(muscInd:end),'Kse'),1,'first')+muscInd-1;
        kpInd = find(contains(inText(muscInd:end),'Kpe'),1,'first')+muscInd-1;
        factor = muscDat{ii,3}/muscDat{ii,2}; 
        if ~(factor > .99 && factor < 1.01)
            origKs = double(extractBetween(string(inText{ksInd}),'>','</'));
            origKp = double(extractBetween(string(inText{kpInd}),'>','</'));
            inText{ksInd} = replaceBetween(inText{ksInd},'>','</',num2str(origKs*factor));
            inText{kpInd} = replaceBetween(inText{kpInd},'>','</',num2str(origKp*factor));
        end
    end
    
    % Turn off the constant joint stimuli
    cStimNames = {'ConstantHip';'ConstantKnee';'ConstantxAnkle'};
    for ii = 1:length(cStimNames)
        cStimInd = find(contains(inText,cStimNames{ii}));
        inText{cStimInd+2} = '<Enabled>False</Enabled>';
    end
    
    testSim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_neutral_test.asim";
    fileID = fopen(testSim,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);

    % NW neutral is [103.58, 116.73, 138.11] in NW coords
    % In [5.1427, 14.504, 21.863] in Animatlab coords [0.089757, 0.25314, 0.38158]
    testObj = design_synergy(testSim);
sdata = processSimData(testSim);
testJM(:,1) = sdata(6).data(:,contains(sdata(6).data_cols,'Hip'));
testJM(:,1) = sdata(6).data(:,contains(sdata(6).data_cols,'Knee'));
testJM(:,1) = sdata(6).data(:,contains(sdata(6).data_cols,'Ankle'));

figure;
subplot(2,1,1)
plot(testObj.theta_motion.*(180/pi)+[98.4373 102.226 116.2473],'LineWidth',3);
subplot(2,1,2)
plot(nobj.theta_motion.*(180/pi)+[98.4373 102.226 116.2473],'LineWidth',1.5);