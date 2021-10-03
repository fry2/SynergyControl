muscleSim = [pwd,'\Animatlab\SynergyWalking\muscleStim.asim'];
% motorProj = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109.aproj';
% motorSim = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim';
motorProj = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced.aproj';
motorSim = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim';

jointAngles2Synergies
% synergyProjectBuilder
[nsys] = indivProjectBuilder(motorProj,motorSim,current2inject,forces,obj);

%motorsim = [pwd,'\Animatlab\SynergyWalking\motorStim.asim'];

%motorsim = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim';
%motorsim = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_motorsim.asim';

muscdata = processSimData(muscleSim);
motordata = processSimData(motorSim);
close all
clear pts mls
numMuscles = length(obj.musc_obj);
    pts = zeros(numMuscles,length(obj.musc_obj{1}.passive_tension));
    for ii = 1:numMuscles
        muscle = obj.musc_obj{ii};
        pts(ii,:) = muscle.passive_tension;
        mls(ii,:) = muscle.muscle_length_profile(obj.sampling_vector,1)';
    end

backer = 10;
plotCell = {'PassiveTension',muscdata(2).data(1:end-backer,:),pts';...
            'JointMotion',muscdata(1).data(1:end-backer,:),motordata(2).data;...
            'MN Stimulation',muscdata(3).data(1:end-backer,:),V_musc';...
            'Muscle Length',muscdata(4).data(1:end-backer,:),mls';...
            'MuscleTension',muscdata(5).data(1:end-backer,:),forces;...
            'Muscle Am',muscdata(6).data(1:end-backer,:),Am_musc';...
            'Muscle Al',muscdata(7).data(1:end-backer,:),100.*Al_musc_all'};

numSubs = size(plotCell,1);
%screensize = [-1919,121,1920,1004];
screensize = [1,301,1440,824];
if 0
    for jj = 1:numSubs
        % Reorder muscles
        if contains(plotCell{jj,1},'Musc')
            muscSimNames = muscdata(jj).data_cols;
            oldSimData = plotCell{jj,2};
            newSimData = zeros(size(oldSimData));
            for ii = 1:numMuscles
                muscObjName = obj.musc_obj{ii}.muscle_name;
                simInd = find(contains(muscSimNames,muscObjName(4:end)));
                if isempty(simInd)
                    keyboard
                end
                newSimData(:,ii) = oldSimData(:,simInd);
            end
            plotCell{jj,2} = newSimData;
        end
        if strcmp(plotCell{jj,1},'Muscle Length')
            meanSim = mean(plotCell{jj,2})';
            meanCalc = mean(mls,2);

            figure('Position',screensize)
            a = subplot(3,1,1);
            bar(meanSim.*1000);
            b = subplot(3,1,2);
            bar(meanCalc.*1000);
            c = subplot(3,1,3);
            bar((meanSim-meanCalc).*1000);
        else
            figure('Position',screensize)
        %     vTemp = plotCell{jj,2};
        %     vTemp = interp1(1:length(vTemp),vTemp,linspace(1,length(vTemp),length(plotCell{jj,3})));
            yMax = max([max(plotCell{jj,2},[],'all') max(plotCell{jj,3},[],'all')]);
            yMin = min([min(plotCell{jj,2},[],'all') min(plotCell{jj,3},[],'all')]);
            subplot(2,1,1)
            plot(plotCell{jj,2},'LineWidth',3)
            ylim([yMin yMax])
            ylabel('Sim Muscle Driven','FontSize',18)
            xlim([0 length(plotCell{jj,2})])
            title(plotCell{jj,1})
            subplot(2,1,2)
            plot(plotCell{jj,3},'LineWidth',3)
            ylabel('Calc Based on Motor','FontSize',18)
            ylim([yMin yMax])
        end
    end
end
        
% jmInd = find(contains({muscdata.name},'JointMotion'));
% bb = [obj.theta_motion ; muscdata(jmInd).data];
% figure
% subplot(2,1,1)
%     plot(obj.theta_motion_time,obj.theta_motion,'LineWidth',3)
%     title('Motor Driven Input Waveforms','FontSize',20)
%     ylabel('Joint Angle','FontSize',1close8)
%     legend({'Hip';'Knee';'Ankle'})
%     ylim([min(bb,[],'all') max(bb,[],'all')])
%         cw = obj.theta_motion(obj.sampling_vector(1):obj.sampling_vector(end),:);
%         sw = [muscdata(jmInd).data(:,3) muscdata(jmInd).data(:,2) muscdata(jmInd).data(:,1)];
%         cw2 = interp1(1:length(cw),cw,linspace(1,length(cw),length(sw)));
% subplot(2,1,2)
%     temp1 = plot(muscdata(2).time,cw2,'LineWidth',3);
%     colorvec = get(temp1,'Color');
%     colorvec = [cell2mat(colorvec),[.2 ;.2; .2]];
%         for ii = 1:3
%             temp1(ii).Color = colorvec(ii,:);
%         end
%     hold on
%     hipInd = find(contains(muscdata(jmInd).data_cols,'Hip'),1,'first');
%     kneeInd = find(contains(muscdata(jmInd).data_cols,'Knee'),1,'first');
%     ankleInd = find(contains(muscdata(jmInd).data_cols,'Ankle'),1,'first');
%     temp2 = plot(muscdata(end).time,[muscdata(jmInd).data(:,hipInd) muscdata(jmInd).data(:,kneeInd) muscdata(jmInd).data(:,ankleInd)],'LineWidth',3);
%         for ii = 1:3
%             temp2(ii).Color = [colorvec(ii,1:3),1];
%         end
%     legend({'Calculated';'Simulated'})
%     title('Muscle Driven Output Waveforms','FontSize',20)
%     ylabel('Joint Angle','FontSize',18)
%     legend({'Hip';'Knee';'Ankle'})
%     ylim([min(bb,[],'all') max(bb,[],'all')])
% 
% obj_fake = design_synergy(sim_file);