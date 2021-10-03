    rootPath = "G:\My Drive\Rat\Swing Paper";
    rootContents = struct2cell(dir(rootPath))';
    rootContents = rootContents(3:end,:);
    dTorque_legender = []; count = 0;

    scaleLookup = {'Rat',1;...
                   'Cat',2.63;...
                   'Horse',12;...
                   'Goat',6.46;...
                   'Human',8.66;...
                   'Cuis',1.1;...
                   'TreeShrew',.9;...
                   'Pika',1;...
                   'StickInsect',.375};
    scaleLookup = sortrows(scaleLookup,2);

    for ii = 1:size(rootContents,1)
        if rootContents{ii,5} == 1 % object is a folder
            dPath = [char(rootPath),'\',char(rootContents{ii,1})];
            dContents = struct2cell(dir(dPath))';
            dContents = dContents(3:end,:);
            dLog = and(contains(dContents(:,1),'JointMoments'),contains(dContents(:,1),'.xlsx'));
            for kk = 1:length(dLog) % Does this folder have a JointMotion.xlsx spreadsheet?
                if dLog(kk) == 1
                    dInd = find(dLog);
                    dFilename = dContents{kk,1}(1:end-5);
                    dInfo = strsplit(dFilename,'_');
                    outData = normalize_excel_data([dContents{kk,2},'\',dContents{kk,1}],dInfo{3});
                    if strcmp(dInfo{4},'Fpos') % My model is Flexion negative. Switch data to match
                        outData = -1.*outData;
                    end
                    slRow = contains(scaleLookup(:,1),dInfo{2});
                    if size(scaleLookup,2) == 2 % No data has been added yet, will throw error
                        scaleLookup{slRow,3} = outData;
                        scaleLookup{slRow,4}{1} = [dInfo{2},' (',rootContents{ii,1},')'];
                    else
                        if isempty(scaleLookup{slRow,3})
                            scaleLookup{slRow,3} = outData;
                            scaleLookup{slRow,4}{1} = [dInfo{2},' (',rootContents{ii,1},')'];
                        else
                            scaleLookup{slRow,3}(:,:,size(scaleLookup{slRow,3},3)+1) = outData;
                            scaleLookup{slRow,4}{size(scaleLookup{slRow,3},3)} = [dInfo{2},' (',rootContents{ii,1},')'];
                        end
                    end
                    count = count + 1;
                end
            end
        end
    end
    dTorque_legender = cell(1,count); dTorque_scaler = zeros(1,count);torqueData = zeros(100,3,count);
    count = 1;
    for ii = 1:size(scaleLookup,1)
        if ~isempty(scaleLookup{ii,3})
            for kk = 1:size(scaleLookup{ii,3},3)
                for jj = 1:3
                    torqueData(1:100,jj,count) = scaleLookup{ii,3}(:,jj,kk);
                end
                dTorque_legender{count} = scaleLookup{ii,4}{kk};
                dTorque_scaler(count) = scaleLookup{ii,2};
                count = count + 1;
            end
        end
    end
clear ii jj kk count dInd dInfo dLog slRow
torqueNorm = zeros(100,3,size(torqueData,3));
for ii = 1:size(torqueData,3)
    if 0
        swingNorm  = -1 + 2.*(torqueData(1:37,:,ii) - min(torqueData(1:37,:,ii)))./(max(torqueData(1:37,:,ii)) - min(torqueData(1:37,:,ii)));
        stanceNorm = -1 + 2.*(torqueData(38:end,:,ii) - min(torqueData(38:end,:,ii)))./(max(torqueData(38:end,:,ii)) - min(torqueData(38:end,:,ii)));
        torqueNorm(:,:,ii) = [swingNorm;stanceNorm];
    else
        swingNorm = torqueData(1:37,:,ii)./max(abs(torqueData(1:37,:,ii)));
        stanceNorm = torqueData(38:end,:,ii)./max(abs(torqueData(38:end,:,ii)));
        torqueNorm(:,:,ii) = [swingNorm;stanceNorm];
        %torqueNorm(:,:,ii) = torqueData(:,:,ii)./max(abs(torqueData(:,:,ii)));
    end
end
%% Plot Processed Animal Data
figure('Position',[961,1,960,1003]);jointNames = {'Hip';'Knee';'Ankle'};
for ii = 1:3
    subplot(3,1,ii)
    plot(squeeze(torqueNorm(:,ii,:)),'LineWidth',2);legend(dTorque_legender,'Location','eastoutside');xline(37,'HandleVisibility','off') % Rat
    xlabel('Stride (%)'); ylabel('Normalized Joint Torque'); title([jointNames{ii},' Joint Torque'],'FontSize',16)
end
%% Plot bar graphs of SWING
stepInds = [1,37,100]; range2plot = stepInds(1):stepInds(3); ft2 = []; showStance = 0;
jointNames = {'Hip';'Knee';'Ankle'};

    if showStance
        phaseInds = [floor((1/3)*(stepInds(3)-stepInds(2)+1)),...
                     floor((2/3)*(stepInds(3)-stepInds(2)+1)),...
                     (stepInds(3)-stepInds(2))];
        strideName = 'Stance';
    else
        phaseInds = [floor((1/3)*(stepInds(2)-stepInds(1))),...
                     floor((2/3)*(stepInds(2)-stepInds(1))),...
                     (stepInds(2)-stepInds(1))];
        strideName = 'Swing';
    end

for joint = 1
    if showStance
        ft2 = torqueNorm((stepInds(2)+1):stepInds(3),:,joint);
    else
        ft2 = torqueNorm(stepInds(1):stepInds(2),:,joint);
    end
    sw2(:,1) = mean(ft2(1:phaseInds(1),:));
    sw2(:,2) = mean(ft2(phaseInds(1):phaseInds(2),:));
    sw2(:,3) = mean(ft2(phaseInds(2):phaseInds(3),:));
end

lens2plot = 1:11;

% Establish a color scale based on a 25 unit spread
%     lengthVals = exp(linspace(log(.05),log(50),24)); 
%     ind = find((1-lengthVals)>0,1,'last');
%     lengthVals = [lengthVals(1:ind),1,lengthVals(ind+1:end)]; numLens = length(lengthVals);
    cm = hsv(length(lengthVals));

figure('Position',[961,1,960,1003],'Name',[jointNames{joint},' Normalized Torque']);
    temp = sw2(lens2plot,:)'; dLog = sum(isnan(temp))~=3; temp = temp(:,dLog); dVals = find(dLog);
    bPlot = bar(temp,'FaceAlpha',1,'EdgeAlpha',0,'FaceColor','flat');hold on;
    % Set scale colors
    for ii = 1:size(temp,2)
       dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{dVals(ii)}),'',' ('))),2};
       [~,minVal] = min(abs(dataScale-lengthVals));
       bPlot(ii).CData = cm(minVal,:); 
    end
    % Set tick data
    set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
    set(gca,'YTick',-1:.5:1,'YTickLabel',{'Flexion  -1';'-.5';'0';'.5';'Extension  1'})
    % Set labels
    legend(dTorque_legender(dLog),'Location','eastoutside'); ylim([-1, 1]); grid on
    ylabel('Normalized Joint Torque'); xlabel('Stride Phase')
    title([jointNames{joint},' Normalized Torque'])
%% function import data from all sheets
function outData = process_excel_data(inSimPath)
    [~,sheets] = xlsfinfo(inSimPath);
    outData = cell(length(sheets),1);
    for ii = 1:length(sheets)
        outData{ii,1} = sheets{ii};
        outData{ii,2} = table2array(readtable(inSimPath,'Sheet',sheets{ii}));
    end
end
%% function: normalize data
function outData = normalize_excel_data(inSimPath,dType)
    rawData = process_excel_data(inSimPath); outData = zeros(100,size(rawData,1));
    fullStride = or(contains(dType,'SwSt'),contains(dType,'StSw'));
    if fullStride
        strideSplit = str2double(dType(5:end));
        if isempty(strideSplit)
            error('Include the stridesplit after "SwSt" or "StSw" in the filename.')
        end
        dType = dType(1:4);
    end
    for ii = 1:size(rawData,1) % For each joint
        if isempty(rawData{ii,2})
            outData(:,ii) = NaN(100,1);
        else
            percent = rawData{ii,2}(:,1); value = rawData{ii,2}(:,2);
            switch dType
                case 'Sw'
                    swing = interp1(1:length(value),value,linspace(1,length(value),37));
                    stance = NaN(1,63);
                case 'St'
                    swing = NaN(1,37);
                    stance = interp1(1:length(value),value,linspace(1,length(value),63));
                case 'SwSt'
                    temp = interp1(percent,value,linspace(1,max(percent),100));
                    swing = interp1(1:strideSplit,temp(1:strideSplit),linspace(1,strideSplit,37));
                    stance = interp1(strideSplit+1:length(temp),temp(strideSplit+1:end),linspace(strideSplit+1,100,63));
                case 'StSw'
                    temp = interp1(percent,value,linspace(1,max(percent),100));
                    stance = interp1(1:strideSplit,temp(1:strideSplit),linspace(1,strideSplit,63));
                    swing = interp1(strideSplit+1:length(temp),temp(strideSplit+1:end),linspace(strideSplit+1,100,37));
            end
            outData(:,ii) = [swing';stance'];
        end
    end
end