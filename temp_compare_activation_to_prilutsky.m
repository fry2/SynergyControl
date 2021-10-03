% Compare muscle activation patterns to those found in Prilutsky 2015 Computing Motion Dependent... Fig. 10.2.b
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_walking.asim";
%inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_catWalking.asim";
outPath = legScaler(inSimPath,2); obj = design_synergy(outPath);
%%
results_cell = obj.pedotti_optimization;
[Am_musc,V_musc,Am_musc_raw] = Am_generator(obj,results_cell{2,2});
%[total_joint_torques,passive_muscle_torque,inertial_torque] = compute_total_joint_torque(obj,0);
%%
pril = [2,33,12,11,11,18,16,28,17,29,22,23];
prilNames = {'IP';'BFA';'RF';'BFP';[];'VI';'VL';'VM';'MG';'LG';'TA';'SO'};
figure('Position',[962,2,958,994]);
stepSeq = [654,1034,1676]; % Toe Off, Touch Down, Toe Off RAT
%stepSeq = [385,724,1056]; % Toe Off, Touch Down, Toe Off CAT
for ii = 1:length(pril)
    if ii ~= 5
        swing = Am_musc(pril(ii),stepSeq(1):stepSeq(2));
        swing2 = interp1(linspace(1,length(swing),length(swing)),swing,linspace(1,length(swing),33));
        stance = Am_musc(pril(ii),stepSeq(2):stepSeq(3));
        stance2 = interp1(linspace(1,length(stance),length(stance)),stance,linspace(1,length(stance),67));
        subplot(length(pril),1,ii)
        plot(1:100,[swing2,stance2],'LineWidth',3)
        %plot(Am_musc(pril(ii),:),'LineWidth',3)
        ylabel(prilNames{ii},'Rotation',0,'FontSize',14,'HorizontalAlignment','right','VerticalAlignment','middle')
        xline(33)
        xlim([0,100]);
    end
end
%% Plot all muscles from a specific group
counter = 1;
figure('Position',[962,2,958,994]);
mTemp = zoning_sorter(inSimPath,6);
groupNum = 5;
numMuscs = sum(cell2mat(mTemp(:,2)) == groupNum);
for ii = 1:length(mTemp)
    if mTemp{ii,2} == groupNum
        subplot(numMuscs,1,counter)
        %plot(Am_musc(ii,stepSeq(1):stepSeq(3))','LineWidth',3)
        
        swing = Am_musc(ii,stepSeq(1):stepSeq(2));
        swing2 = interp1(linspace(1,length(swing),length(swing)),swing,linspace(1,length(swing),33));
        stance = Am_musc(ii,stepSeq(2):stepSeq(3));
        stance2 = interp1(linspace(1,length(stance),length(stance)),stance,linspace(1,length(stance),67));
        plot(1:100,[swing2,stance2],'LineWidth',3)
        ylabel(obj.musc_obj{ii}.muscle_name(4:end))
%         xlim([0 length(Am_musc(ii,stepSeq(1):stepSeq(3))')])
        xlim([0 100])
        xline(33)
        counter  = counter + 1;
    end
end