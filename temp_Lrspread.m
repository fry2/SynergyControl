% The goal of this script is to determine the effect that a "spread" has on the resting length of muscles
% In this context, a spread is a precent change from a neutral length, as determined by setting the leg to a neutral position
% Originally, resting lengths were set equal to whatever the muscle length is when the elg is in a neutral position
% The problem this generates is that entire groups of muscles switch on and off as the leg passes through this neutral position
% This leads sharp discontinuities in the passive torque of each joint
% To combat this, we implement a random spread (+/- some percent) from the neutral length for eahc muscle
% This script cycles through a number of different spreads to determine how the resulting passive torque profile looks

neutralPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_neutral.asim";
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim";
trialSavePath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_lrspread.asim";
rVals = linspace(.05,.5,10); pTorque = cell(1,length(rVals));
count = 1;
%% Change LT values
    for jj = 1:length(rVals)
        neutral_lengths = LRsolver(neutralPath,rVals(jj));
        inText = importdata(inSimPath);
        mInds = find(contains(inText,'<Type>LinearHillMuscle</Type>'))-2;
        for ii = 1:length(mInds)
            mInd = mInds(ii);
            mName = lower(char(extractBetween(string(inText{mInd}),'>','<')));
            dataInd = find(contains(neutral_lengths.data(:,1),mName(4:end)));
            dataVec = neutral_lengths.data(dataInd,:);
            ltInd = find(contains(inText(mInd:end),'<LengthTension>'),1,'first')+mInd-1;
            inText{ltInd+4} = replaceBetween(inText{ltInd+4},'>','<',num2str(dataVec{4}));
            inText{ltInd+5} = replaceBetween(inText{ltInd+5},'>','<',num2str(dataVec{5}));
            inText{ltInd+8} = replaceBetween(inText{ltInd+8},'>','<',num2str(dataVec{2}));
            inText{ltInd+9} = replaceBetween(inText{ltInd+9},'>','<',num2str(dataVec{3}));
        end

        fileID = fopen(trialSavePath,'w');
        fprintf(fileID,'%s\n',inText{:});
        fclose(fileID);

        [temph,hjm,tempha,moH,mlH,hObj] = muscle_to_joint_VE(trialSavePath,1);
        [tempk,kjm,tempka,moK,mlK,kObj] = muscle_to_joint_VE(trialSavePath,2);
        [tempa,ajm,tempaa,moA,mlA,aObj] = muscle_to_joint_VE(trialSavePath,3);

        pTorqueh = compute_passive_joint_torque(hObj);
        pTorquek = compute_passive_joint_torque(kObj);
        pTorquea = compute_passive_joint_torque(aObj);
        pTorque{1,jj} = [-pTorqueh(1,:)',-pTorquek(2,:)',-pTorquea(3,:)'];
        disp(['Finished ',num2str(count),'/',num2str(length(rVals))]);
        count = count + 1;
    end
%%
figure('Position',[962,2,958,994]);
rvinds = 1:3:length(rVals);
legVec = cellfun(@num2str,num2cell(rVals(rvinds)),'UniformOutput',false); cm = winter(length(rVals));

jointMotion = [hObj.theta_motion(100:end-100,1),kObj.theta_motion(100:end-100,2),aObj.theta_motion(100:end-100,3)].*(180/pi)+[98.4373 102.226 116.2473];
s1 = subplot(4,1,1);
plot(hObj.theta_motion_time(100:end-100),jointMotion,'LineWidth',3)
ylabel('Joint Angle (deg)')
xlabel('Time (s)')
title('Effect of Varying Resting Lengths About the Neutral Position','FontSize',16)
for jj = 1:length(rvinds)
    s2 = subplot(4,1,2);
        plot(hObj.theta_motion_time(100:end-100),pTorque{rvinds(jj)}(100:end-100,1),'Color',cm(rvinds(jj),:),'LineWidth',3); hold on
        xlabel('Time (s)')
        ylabel('Torque (Nm)')
    s3 = subplot(4,1,3);
        plot(kObj.theta_motion_time(100:end-100),pTorque{rvinds(jj)}(100:end-100,2),'Color',cm(rvinds(jj),:),'LineWidth',3); hold on
        xlabel('Joint Angle (RAD)')
        ylabel('Torque (Nm)')
    s4 = subplot(4,1,4);
        plot(aObj.theta_motion_time(100:end-100),pTorque{rvinds(jj)}(100:end-100,3),'Color',cm(rvinds(jj),:),'LineWidth',3); hold on
        xlabel('Joint Angle (RAD)')
        ylabel('Torque (Nm)')
end
    legend(s1,{'Hip';'Knee';'Ankle'},'Location','eastoutside')
    legend(s2,legVec,'Location','eastoutside')
    legend(s3,legVec,'Location','eastoutside')
    legend(s4,legVec,'Location','eastoutside')