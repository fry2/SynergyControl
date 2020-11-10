%%
projPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl.aproj';
simPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim';

%%
% simfile = importdata(simPath);
% kneebounds = [-.6 2*.05155];
% kneecoeffs = [sum(kneebounds)/2 sum(abs(kneebounds))/2];
% %kneecoeffs_signs = sign(kneecoeffs);
% %kneecoeffs = cellfun(@num2str,num2cell(abs(kneecoeffs)),'UniformOutput',false);
% eqInd = find(contains(simfile,'<Name>zMaxMinKnee</Name>'))+10;
% newKneeEq = ['<Equation>0,',num2str(abs(kneecoeffs(1))),',',num2str(kneecoeffs(2)),',t,sin,*,+,-</Equation>'];
% simfile{eqInd} = newKneeEq;
%         fileID = fopen(simPath,'w');
%         fprintf(fileID,'%s\n',simfile{:});
%         fclose(fileID);
%%
obj = design_synergy(simPath);
results_cell = pedotti_optimization(obj);
%%
fa = results_cell{3,2}';
[nsys] = indivProjectBuilder(projPath,simPath,results_cell{2,2}',obj);
sim_file_revised = strcat(simPath(1:end-5),'_fake.asim');
sdata = processSimData(sim_file_revised);
clear force_jm
for ii = 1:3
    switch ii
        case 1
            joint = 'Hip';
        case 2
            joint = 'Knee';
        case 3
            joint = 'Ankle';
    end
    force_jm(:,ii) = sdata(1).data(:,find(contains(sdata(1).data_cols,joint)));
end
force_jm = force_jm.*(180/pi)+[98.4373 102.226 116.2473];
%%
ef = results_cell{9,2};
t1 = find(gradient(ef)==1.5);
t1 = t1(1:2:end);
t2 = find(gradient(ef)==-1.5);
t2 = t2(1:2:end);
xliner = sort([t1,t2]);
figure('Position',[962,2,958,994]);
minlen = min([length(obj.theta_motion),length(force_jm)]);
shifted_jm = obj.theta_motion(1:minlen,:).*(180/pi)+[98.4373 102.226 116.2473];
force_jm = force_jm(1:minlen,:);
jointLims = [min([shifted_jm;force_jm],[],'all'), max([shifted_jm;force_jm],[],'all')];
actLims = 1.1.*[min([results_cell{2,2}(:,335:end)';sdata(5).data(335:end,:)],[],'all'), max([results_cell{2,2}(:,335:end)';sdata(5).data(335:end,:)],[],'all')];
time = obj.theta_motion_time(1:minlen);
titlesize = 15; axsize = 12; legsize = 12;
subplot(4,1,1)
    plot(time,shifted_jm,'LineWidth',3)
    for ii = 1:length(xliner)
        xline(xliner(ii))
    end
    title('Desired Joint Motion','FontSize',titlesize)
    ylim(jointLims)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Joint Angle (deg)','FontSize',axsize)
    xlim([0 max(time)])
    legend({'Hip';'Knee';'Ankle'},'Location','southwest','FontSize',legsize)
subplot(4,1,2)
    plot(time,force_jm,'LineWidth',3)
    title('Force-Driven Simulation Joint Motion','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Joint Angle (deg)','FontSize',axsize)
    xlim([0 max(time)])
    ylim(jointLims)
    legend({'Hip';'Knee';'Ankle'},'Location','southwest','FontSize',legsize)
subplot(4,1,3)
    plot(time,results_cell{2,2}(:,1:minlen)','LineWidth',1.5)
    title('Desired Forces','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Muscle Force (N)','FontSize',axsize)
    ylim(actLims)
    xlim([0 max(time)])
subplot(4,1,4)
    plot(time,sdata(5).data(1:minlen,:),'LineWidth',1.5)
    title('Force-Driven Simulation Forces','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Muscle Force (N)','FontSize',axsize)
    ylim(actLims)
    xlim([0 max(time)])
    return
%% Plot the muscle forces and rank the "worst offender" muscles 
jointDiff = sdata(5).data(1:length(results_cell{2,2}),:)-results_cell{2,2}';
figure;plot(jointDiff);
for ii = 1:38
    bigdiffs{ii,1} = obj.musc_obj{ii}.muscle_name;
    bigdiffs{ii,2} = jointDiff(5765,ii);
end
bigdiffs = sortrows(bigdiffs,2,'descend');