function [newTens,origTens,newJM,origJM] = paramChanger(inSimPath,scales,sdata)
    % For a given simulation file, scale all muscle muscle parameters by a certain value
    % inSimPath: char: simulation path
    % scales: double (1x3): the amount to scale [B, Ks, Kp] % ex: [1 1 10]
    %sdata = processSimData(inSimPath);
    tensionInd = contains({sdata.name},'KeyMuscTen');
    origTens = sdata(tensionInd).data;
    
    jmInd = contains({sdata.name},'JointMotion');
    jointNames = {'Hip';'Knee';'Ankle'}; jointOffsets = [98.4373 102.226 116.2473]; origJM = zeros(length(sdata(jmInd).data),3);
    for ii = 1:3
        jInd = find(contains(sdata(jmInd).data_cols,jointNames{ii}));
        origJM(:,ii) = sdata(jmInd).data(:,jInd).*(180/pi)+jointOffsets(ii);
    end
    
    origtext = importdata(inSimPath);
    muscInds = find(contains(origtext,'<Type>LinearHillMuscle</Type>'));
    for ii = 1:length(muscInds)
        % B
        bInd = find(contains(origtext(muscInds(ii):end),'<B>'),2)+muscInds(ii)-1;
        bInd = bInd(2);
        oldPar = double(extractBetween(string(origtext{bInd}),'>','<'));
        origtext{bInd} = replaceBetween(origtext{bInd},'>','<',num2str(oldPar*scales(1)));
        % Ks
        ksInd = find(contains(origtext(muscInds(ii):end),'<Kse>'),1)+muscInds(ii)-1;
        oldPar = double(extractBetween(string(origtext{ksInd}),'>','<'));
        origtext{ksInd} = replaceBetween(origtext{ksInd},'>','<',num2str(oldPar*scales(2)));
        % Kp
        kpInd = find(contains(origtext(muscInds(ii):end),'<Kpe>'),1)+muscInds(ii)-1;
        oldPar = double(extractBetween(string(origtext{kpInd}),'>','<'));
        origtext{kpInd} = replaceBetween(origtext{kpInd},'>','<',num2str(oldPar*scales(3)));
        % Fmax
        fInd = find(contains(origtext(muscInds(ii):end),'<MaximumTension>'),1)+muscInds(ii)-1;
        oldPar = double(extractBetween(string(origtext{fInd}),'>','<'));
        origtext{fInd} = replaceBetween(origtext{fInd},'>','<',num2str(oldPar*scales(4)));
    end
    
    txtInd1 = find(contains(origtext,'OutputFilename>Joint'));
    origtext{txtInd1} = replaceBetween(origtext{txtInd1},'ion','.txt','_temp');
    txtInd2 = find(contains(origtext,'OutputFilename>KeyMuscTen'));
    origtext{txtInd2} = replaceBetween(origtext{txtInd2},'Ten','.txt','_temp');
    
    jobSavePath = [pwd,'\temp1.asim'];
    fileID = fopen(jobSavePath,'w');
    fprintf(fileID,'%s\n',origtext{:});
    fclose(fileID);
    
    newData = processSimData(jobSavePath);
    tensionInd2 = contains({newData.name},'KeyMuscTen');
    newTens = newData(tensionInd2).data;
    
    jmInd = contains({newData.name},'JointMotion');
    newJM = zeros(length(newData(jmInd).data),3);
    for ii = 1:3
        jInd = find(contains(newData(jmInd).data_cols,jointNames{ii}));
        newJM(:,ii) = newData(jmInd).data(:,jInd).*(180/pi)+jointOffsets(ii);
    end
   
    muscZones = zoning_sorter(inSimPath,6);
    cm = [0,1,0;1,0,0;0.9216,0.8235,0.2039;0,0,1;1,0,1;0,1,1];
    
    figure('Position',[962,2,958,994]);
    titlesize = 15; axsize = 12; legsize = 12;
    subplot(4,1,1)
        plot(origJM,'LineWidth',3)
        title('Original Joint Motion','FontSize',titlesize)
        xlabel('Time (s)','FontSize',axsize)
        ylabel('Joint Angle (deg)','FontSize',axsize)
        legend({'Hip';'Knee';'Ankle'},'Location','southeast','FontSize',legsize)
    subplot(4,1,2)
        plot(newJM,'LineWidth',3)
        title(['Joint Motion w/ Scaled Values ','B= ',num2str(scales(1)),...
            ', Ks= ',num2str(scales(2)),', Kp= ',num2str(scales(3)),', Fmax= ',num2str(scales(4))],'FontSize',titlesize)
        xlabel('Time (s)','FontSize',axsize)
        ylabel('Joint Angle (deg)','FontSize',axsize)
        legend({'Hip';'Knee';'Ankle'},'Location','southeast','FontSize',legsize)
    subplot(4,1,3)
        for ii = 1:size(origTens,2)
            plot(origTens(:,ii),'LineWidth',1.5,'Color',cm(muscZones{ii,2},:))
            hold on
        end
        title('Original Forces','FontSize',titlesize)
        xlabel('Time (s)','FontSize',axsize)
        ylabel('Muscle Force (N)','FontSize',axsize)
    subplot(4,1,4)
        for ii = 1:size(newTens,2)
            plot(newTens(:,ii),'LineWidth',1.5,'Color',cm(muscZones{ii,2},:))
            hold on
        end
        title('Forces w/ Scaled Values','FontSize',titlesize)
        xlabel('Time (s)','FontSize',axsize)
        ylabel('Muscle Force (N)','FontSize',axsize)
        
%     % Figure with each muscle zone in its own subplot
%     figure('Position',[962,2,958,994]);
%     titlesize = 15; axsize = 12; legsize = 12; linew = 1;
%     zoneTens = cell(6);
%     for ii = 1:38
%         zoneTens{muscZones{ii,2}}(:,end+1) = newTens(:,ii);
%     end
%     subplot(6,1,1)
%         plot(zoneTens{1},'Color',cm(1,:),'LineWidth',linew)
%     subplot(6,1,2)
%         plot(zoneTens{2},'Color',cm(2,:),'LineWidth',linew)
%     subplot(6,1,3)
%         plot(zoneTens{3},'Color',cm(3,:),'LineWidth',linew)
%     subplot(6,1,4)
%         plot(zoneTens{4},'Color',cm(4,:),'LineWidth',linew)
%     subplot(6,1,5)
%         plot(zoneTens{5},'Color',cm(5,:),'LineWidth',linew)
%     subplot(6,1,6)
%         plot(zoneTens{6},'Color',cm(6,:),'LineWidth',linew)

% Put this into command line to generate images for a range of scales for a parameter
%  scaleVals = linspace(.05,5,10);
% for ii = 1:length(scaleVals)
% [newTens,origTens,newJM,origJM] = paramChanger(sim_file_revised,[1,scaleVals(ii),1,1],sdata);
% saveFolder = 'G:\My Drive\Rat\MeetingFiles\Meeting_20210125\KsScale\';
% integ=floor(scaleVals(ii));
% fract=scaleVals(ii)-integ;
% saveas(gcf,[saveFolder,'Ks',num2str(integ),'p',num2str(round(fract,3)),'.png'])
% close(gcf)
% end
    
    delete(jobSavePath)
end