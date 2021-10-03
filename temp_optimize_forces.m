%%
projPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl.aproj';
simPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim';

%%
%     simfile = importdata(simPath);
%     fmaxinds = find(contains(simfile,'<MaximumTension>'));
%     for ii = 1:length(fmaxinds)
%         oldNum = double(extractBetween(string(simfile{fmaxinds(ii)}),'>','</'));
%         simfile{fmaxinds(ii)} = ['<MaximumTension>',num2str(10*oldNum),'</MaximumTension>'];
%     end
%     fileID = fopen(simPath,'w');
%     fprintf(fileID,'%s\n',simfile{:});
%     fclose(fileID);
%% Create the FullLeg object and run the force optimization
obj = design_synergy(simPath);
results_cell = pedotti_optimization(obj);
%% Generate a force-driven project
%fa = results_cell{3,2};
nsys = indivProjectBuilder(projPath,simPath,results_cell{2,2},obj);
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

figure('Position',[962,2,958,994]);
minlen = min([length(obj.theta_motion),length(force_jm)]);
time = obj.theta_motion_time(1:minlen);
% Sort out the joint angle waveforms
    shifted_jm = obj.theta_motion(1:minlen,:).*(180/pi)+[98.4373 102.226 116.2473];
    force_jm = force_jm(1:minlen,:);
    jointLims = [min([shifted_jm;force_jm],[],'all'), max([shifted_jm;force_jm],[],'all')];
% Sort out the individual force waveforms
    desiredForces = results_cell{2,2}(1:minlen,:);
    simForces = sdata(5).data(1:minlen,:);
    actLims = 1.05.*[min([desiredForces(1:end-300,:);simForces(1:end-300,:)],[],'all'), max([desiredForces(1:end-300,:);simForces(1:end-300,:)],[],'all')];
% Sort out error shading
    ef = results_cell{9,2};
    has_errors = sum(ef) ~= length(ef);
    if has_errors
        maxJM = max(shifted_jm,[],'all'); minJM = min(shifted_jm,[],'all'); maxEF = max(ef); minEF = min(ef);
        amp = (maxJM-minJM)/(maxEF-minEF); shift = (maxJM+minJM)/2-((maxEF+minEF)/2)*amp;
        x = time'; y1 = maxJM.*ones(minlen,1)'; y2 = (amp.*ef(1:minlen)+shift)';
    end
% Sort out force coloring
    muscZones = zoning_sorter(simPath,6);
    cm = [0,1,0;1,0,0;0.9216,0.8235,0.2039;0,0,1;1,0,1;0,1,1];
  
titlesize = 15; axsize = 12; legsize = 12;
subplot(4,1,1)
    plot(time,shifted_jm,'LineWidth',3)
    hold on
    if has_errors
        patch([x fliplr(x)], [y1 fliplr(y2)], 'r','FaceAlpha',.2,'EdgeAlpha',0)
    end
    title('Desired Joint Motion','FontSize',titlesize)
    %ylim(jointLims)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Joint Angle (deg)','FontSize',axsize)
    xlim([0 max(time)])
    legend({'Hip';'Knee';'Ankle'},'Location','southeast','FontSize',legsize)
subplot(4,1,2)
    plot(time,force_jm,'LineWidth',3)
    title('Force-Driven Simulation Joint Motion','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Joint Angle (deg)','FontSize',axsize)
    xlim([0 max(time)])
    %ylim(jointLims)
    legend({'Hip';'Knee';'Ankle'},'Location','southeast','FontSize',legsize)
subplot(4,1,3)
    for ii = 1:size(desiredForces,2)
        plot(time,desiredForces(:,ii),'LineWidth',1.5,'Color',cm(muscZones{ii,2},:))
        hold on
    end
    title('Desired Forces','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Muscle Force (N)','FontSize',axsize)
    ylim(actLims)
    xlim([0 max(time)])
subplot(4,1,4)
    for ii = 1:size(desiredForces,2)
        plot(time,simForces(:,ii),'LineWidth',1.5,'Color',cm(muscZones{ii,2},:))
        hold on
    end
    title('Force-Driven Simulation Forces','FontSize',titlesize)
    xlabel('Time (s)','FontSize',axsize)
    ylabel('Muscle Force (N)','FontSize',axsize)
    ylim(actLims)
    xlim([0 max(time)])
    return
%% Plot the muscle forces and rank the "worst offender" muscles
minLen = min([length(results_cell{2,2}), length(sdata(5).data)]);
jointDiff = sdata(5).data(1:minLen,:)-results_cell{2,2}(1:minLen,:);
muscZones = zoning_sorter(simPath,6);
%figure;plot(jointDiff);
for mnum = 1:38
    musc = obj.musc_obj{mnum};
    b = musc.damping; ks = musc.Kse; kp = musc.Kpe;
    bigdiffs{mnum,1} = obj.musc_obj{mnum}.muscle_name;
    bigdiffs{mnum,2} = muscZones{mnum,2};
    bigdiffs{mnum,3} = sum(jointDiff(:,mnum).^2);
    bigdiffs{mnum,4} = obj.musc_obj{mnum}.damping/(obj.musc_obj{mnum}.Kse+obj.musc_obj{mnum}.Kpe);% muscle time constant in (s). big == slow
    bigdiffs{mnum,5} = ks/kp;
end
bigdiffs = sortrows(bigdiffs,3,'descend');
%%  Plot different muscle values vs force error. Includes zone coloring
    figure
    cm = [0,1,0;...
          1,0,0;...
          0.9216,0.8235,0.2039;...
          0,0,1;...
          1,0,1;...
          0,1,1];
    for ii = 1:38
        xval = bigdiffs{ii,5};
        scatter(xval,bigdiffs{ii,3},36,cm(bigdiffs{ii,2},:),'o','LineWidth',10)
        hold on
    end
%% 
desiredTdot = smoothdata(gradient(results_cell{2,2},obj.dt_motion),'gaussian',150);
simulatedTdot = smoothdata(gradient(sdata(5).data(1:minlen,:),obj.dt_motion),'gaussian',150);
figure;
subplot(2,1,1)
plot(desiredTdot)
subplot(2,1,2)
plot(simulatedTdot)