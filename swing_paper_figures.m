%% Model Info    
ft_model_tot = []; ft_data = []; legender = cell(1,1); showStance = 1;
  
% Define phase indices
    if showStance
        strideName = 'Stance';
        phaseInds =  [21,42,63];
        range2plot = 38:100;
    else
        strideName = 'Swing';
        phaseInds  = [12,24,37];
        range2plot = 1:37;
    end
%
phase_model = zeros(numLens,3,3);
phase_data = zeros(size(torqueNorm,3),3,3);
for joint = 1:3
        ft_model = torques_act{joint}(range2plot,:);
        ft_data = squeeze(torqueNorm(range2plot,joint,:));
    phase_model(:,:,joint) = [mean(ft_model(1:phaseInds(1),:))',mean(ft_model(phaseInds(1)+1:phaseInds(2),:))',mean(ft_model(phaseInds(2)+1:phaseInds(3),:))'];
    phase_data(:,:,joint)  = [mean(ft_data(1:phaseInds(1),:))',mean(ft_data(phaseInds(1)+1:phaseInds(2),:))',mean(ft_data(phaseInds(2)+1:phaseInds(3),:))'];
end
disp('Model and Data structures made')
%% 3x2 - Model and Data - Bar Graphs Figure
cm = jet(length(lengthVals));
figure('Position',[762,-7,1159,1003],'Name',strideName); jointNames = {'Hip';'Knee';'Ankle'};
    for jj = 1:6
        aa = subplot(3,2,jj);
        jointNum = floor(.5*(jj+1));
        if mod(jj,2)
%             modell2p = [1,6,9,11,14,18,19,24];
            modell2p = 1:numLens;
            %modell2p = [1,6,9,11,14,15,17,18,19];
            bPlot = bar(phase_model(modell2p,:,jointNum)','FaceAlpha',1,'EdgeAlpha',0,'FaceColor','flat');hold on;
            % Set scale colors
            for ii = 1:length(modell2p)
               bPlot(ii).CData = cm(modell2p(ii),:); 
            end
            % Tweak tick info on axes
            set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
            % Add labels and title
            legend(flip(bPlot),flip(cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false)),'Location','eastoutside'); ylim([-1, 1]); grid on
            aa.Position(3:4) = [0.2454,0.2157];
        else
            % Data Info
            numDataLens = 1:size(phase_data,1);
            temp = phase_data(numDataLens,:,jointNum)'; dLog = sum(isnan(temp))~=3; temp = temp(:,dLog); dVals = find(dLog);
            cPlot = bar(temp,'FaceAlpha',1,'EdgeAlpha',0,'FaceColor','flat');hold on;
            % Set scale colors
            for ii = 1:size(temp,2)
               dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{dVals(ii)}),'',' ('))),2};
               [~,minVal] = min(abs(dataScale-lengthVals(modell2p)));
               cPlot(ii).CData = cm(modell2p(minVal),:); 
            end
            legend(flip(cPlot),flip(dTorque_legender(dLog)),'Location','eastoutside'); ylim([-1, 1]); grid on
            aa.Position(3:4) = [0.2454,0.2157];
        end
        set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
        ylabel('Normalized Joint Torque','FontSize',14); xlabel('Stride Phase','FontSize',14)
        title(['Total Joint Torques of the ',jointNames{jointNum}],'FontSize',16)
    end
%% 3x2 - Model and Data - Line Graphs Figure
cm = turbo(length(lengthVals)); jointNames = {'Hip';'Knee';'Ankle'};
figure('Position',[762,-7,1159,1003],'Name',strideName); subSize = [];
    for jj = 1:6
        aa = subplot(3,2,jj);
        jointNum = floor(.5*(jj+1));
        if mod(jj,2)
            modell2p = 1:numLens;
            yline(0,'LineWidth',2,'Alpha',.1); hold on
            bPlot = plot(phase_model(modell2p,:,jointNum)','-s','LineWidth',3);
            % Set scale colors
            for ii = 1:length(modell2p)
               bPlot(ii).Color = cm(modell2p(ii),:); 
               bPlot(ii).MarkerFaceColor = cm(modell2p(ii),:);
            end
            % Tweak tick info on axes
            set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
            % Add labels and title
            legend(flip(bPlot),flip(cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false)),'Location','eastoutside'); ylim([-1, 1]); grid on
            aa.Position(3:4) = [0.2454,0.2157];
        else
            % Data Info
            numDataLens = 1:size(phase_data,1);yline(0,'LineWidth',2,'Alpha',.1); hold on
            temp = phase_data(numDataLens,:,jointNum)'; dLog = sum(isnan(temp))~=3; temp = temp(:,dLog); dVals = find(dLog);
            cPlot = plot(temp,'-s','LineWidth',3);
            % Set scale colors
            for ii = 1:size(temp,2)
               dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{dVals(ii)}),'',' ('))),2};
               [~,minVal] = min(abs(dataScale-lengthVals(modell2p)));
               cPlot(ii).Color = cm(modell2p(minVal),:); 
               cPlot(ii).MarkerFaceColor = cm(modell2p(minVal),:); 
            end
            legend(flip(cPlot),flip(dTorque_legender(dLog)),'Location','eastoutside'); ylim([-1, 1]); grid on
            aa.Position(3:4) = [0.2454,0.2157];
        end
        set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
        ylabel('Normalized Joint Torque','FontSize',14); xlabel('Stride Phase','FontSize',14)
        title(['Total Joint Torques of the ',jointNames{jointNum}],'FontSize',16)
    end
%% 2x3 - Model and Data - Line Graphs Figure
cm = jet(length(lengthVals)); jointNames = {'Hip';'Knee';'Ankle'}; modell2p = 1:numLens;
figure('Position',[1,1,1920,1003],'Name',strideName); subSize = []; dLog = [];
    for jj = 1:6
        aa = subplot(2,3,jj);
%         jointNum = floor(.5*(jj+1));
        yline(0,'LineWidth',3,'Alpha',.15); hold on
        if jj > 3
            % Add labels and title
            if showStance == 1
                if jj ~= 6
                    bPlot = plot(phase_model(modell2p,:,jj-3)','-s','LineWidth',3);
                    % Set scale colors
                    for ii = 1:length(modell2p)
                       bPlot(ii).Color = cm(modell2p(ii),:); 
                       bPlot(ii).MarkerFaceColor = cm(modell2p(ii),:);
                    end
                else
                    bPlot = plot(phase_model(2:6,:,jj-3)','-s','LineWidth',3);
                    % Set scale colors
                    for ii = 1:5
                       bPlot(ii).Color = cm(modell2p(ii+1),:); 
                       bPlot(ii).MarkerFaceColor = cm(modell2p(ii+1),:);
                    end
                end
            else
                bPlot = plot(phase_model(modell2p,:,jj-3)','-s','LineWidth',3);
                % Set scale colors
                for ii = 1:length(modell2p)
                   bPlot(ii).Color = cm(modell2p(ii),:); 
                   bPlot(ii).MarkerFaceColor = cm(modell2p(ii),:);
                end
            end
            if jj == 5
                legend(flip(bPlot),flip(cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false)),'Location','eastoutside',...
                    'Position',[0.82,0.207,0.039,0.15],'FontSize',12);
                text(6.150762463343108,0.602375366568915,1,{'Modeled';'Animal Sizes'},'HorizontalAlignment','center','FontSize',12)
            end
            %aa.Position(3:4) = [0.2454,0.2157];
            title([jointNames{jj-3},' Torque - Modeled'],'FontSize',20)
        else
            % Data Info
            numDataLens = 1:size(phase_data,1);
            temp = phase_data(numDataLens,:,jj)'; dLog(jj,:) = sum(isnan(temp))~=3; temp = temp(:,logical(dLog(jj,:))); dVals = find(dLog(jj,:));
            tempPlot = plot(temp,'-s','LineWidth',3);
            cPlot{jj} = tempPlot;
            clear tempPlot
            % Set scale colors
            for ii = 1:size(temp,2)
               dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{dVals(ii)}),'',' ('))),2};
               [~,minVal] = min(abs(dataScale-lengthVals(modell2p)));
               cPlot{jj}(ii).Color = cm(modell2p(minVal),:); 
               cPlot{jj}(ii).MarkerFaceColor = cm(modell2p(minVal),:); 
            end
            if jj == 3
                plot2leg = find(sum(dLog,2)==max(sum(dLog,2))); % pick a plot to make the legend from. Choose the first one that includes all the lines
                legend(flip(cPlot{plot2leg(1)}),flip(dTorque_legender(any(dLog))),'FontSize',12,'Location','eastoutside','Position',[0.82,0.69,0.084,0.135]);
                text(3.493870967741935,0.713812316715543,1,'Animal Data','FontSize',12)
                clear plot2leg
            end
            %aa.Position(3:4) = [0.2454,0.2157];
            title([jointNames{jj},' Torque - Data'],'FontSize',20)
        end
        ylim([-1 1]); grid on; pbaspect([1 1 1])
        set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
        switch jj
            case 2
                aa.Position = [0.365,0.58,0.21,0.34];
            case 5
                aa.Position = [0.365,0.11,0.21,0.34];
            case 3
                aa.Position = [0.6,0.58,0.21,0.34];
            case 6
                aa.Position = [0.6,0.11,0.21,0.34];
        end
        if jj == 1 || jj == 4
            ylabel('Normalized Joint Torque','FontSize',18);
        elseif jj == 2 || jj == 5
            xlabel('Stride Phase','FontSize',18)
        end
    end
%% 1x3 - Model Only - Line Graphs Figure ALL
cm = jet(length(lengthVals)); jointNames = {'Hip';'Knee';'Ankle'};
    % For defining a custom color map between two colors
    color1 = [0,0,1]; color2 = [1,0,0];
    r = linspace(color1(1),color2(1),numLens)'; g = linspace(color1(2),color2(2),numLens)'; b = linspace(color1(3),color2(3),numLens)';  cm = [r,g,b];
figure('Position',[602,508,1318,481],'Name',strideName); subSize = []; bPlot = [];
    tiledlayout(1,3,'TileSpacing','compact')
    for jj = 1:3
        aa(jj) = subplot(1,3,jj);
            modell2p = 1:numLens; 
            yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on
             bPlot = plot(phase_model(modell2p,:,jj)','-s','LineWidth',3);
            %bPlot = plot(tot_torques(modell2p,:,jj)','-s','LineWidth',3);
            % Set scale colors
            for ii = 1:length(modell2p)
               bPlot(ii).Color = cm(modell2p(ii),:); 
               bPlot(ii).MarkerFaceColor = cm(modell2p(ii),:);
            end
            % Tweak tick info on axes
            set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]},'FontSize',10)
            % Add labels and title
            ylim([-1, 1]); grid on; pbaspect([1,1,1])
            set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
            set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
            switch jj
                case 1
                    ylabel('Normalized Joint Torque','FontSize',16);
                case 2
                    xlabel('Stride Phase','FontSize',16,'Position',[2,-1.675,-1]);
                    title(['Modeled Active Joint Torques - ',strideName],'FontSize',16,'Position',[2,1.2,0])
                case 3
                    lg = legend(flip(bPlot),flip(cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false)),'Location','eastoutside','FontSize',12);
            end
            text(2.6,.9,1,jointNames{jj},'Color','k','FontSize',14);
    end
    
    %     lg = legend(aa(2), cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false),'Location','southoutside',...
    %         'Position',[.3,.028,.4347,.0478],'Orientation','horizontal','FontSize',12);
    
    text(3.6462,.6762,{'Size Relative';'to Rat'},'HorizontalAlignment','center')
    
    for jj = 1:3
        switch jj
            case 1
                set(aa(jj),'Position',[.0716,.1079,.1973,.815])
            case 2
                set(aa(jj),'Position',[.3759,.1079,.1973,.815])
            case 3
                set(aa(jj),'Position',[.6805,.1079,.1973,.815])
        end
    end
    
    miniData = mini_legend_info();
    if showStance
        miniLegSeq = [3,2,1,6,5,4,9,8,7];
    else
        miniLegSeq = 1:9;
    end
    for ii = 1:9
        ind = miniLegSeq(ii);   
        axes('pos',miniData{ii,1}); 
        [YourImage, ~, ImageAlpha] = imread(miniData{ind,2});
        h = imshow(YourImage);
        h.AlphaData = ImageAlpha;
    end
%% 1x3 - Model Only - Line Graphs Figure SIMPLE
cm =     [0 0 1; .5 0 .5; 1 0 0]; jointNames = {'Hip';'Knee';'Ankle'};
figure('Position',[602,508,1318,481],'Name',strideName); subSize = []; bPlot = [];
    tiledlayout(1,3,'TileSpacing','compact')
    for jj = 1:3
        aa(jj) = subplot(1,3,jj);
            switch jj
                case 1 % Hip
                    if showStance == 0 % Is Swing
                        modell2p = [1,3,6];
                    else
                        modell2p = [1,4,6];
                    end
                case 2 % Knee
                    if showStance == 0 % Is Swing
                        modell2p = [1,2,6];
                    else
                        modell2p = [1,2,6];
                    end
                case 3 % Ankle
                    if showStance == 0 % Is Swing
                        modell2p = [1,5,6];
                    else
                        modell2p = [1,3,6];
                    end
            end
            yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on
             bPlot = plot(phase_model(modell2p,:,jj)','-s','LineWidth',3);
            % Set scale colors
            for ii = 1:length(modell2p)
               bPlot(ii).Color = cm(ii,:); 
               bPlot(ii).MarkerFaceColor = cm(ii,:);
            end
            % Tweak tick info on axes
            set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]},'FontSize',10)
            % Add labels and title
            ylim([-1, 1]); grid on; pbaspect([1,1,1])
            set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
            set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
%             lg = legend(bPlot,cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false),...
%                 'Orientation','horizontal','Location','northoutside','FontSize',12);
            switch jj
                case 1
                    ylabel('Normalized Joint Torque','FontSize',16);
                case 2
                    xlabel('Stride Phase','FontSize',16,'Position',[2,-1.675,-1]);
                    title(['Active Muscle Torque During ',strideName,' - Modeled'],'FontSize',16,'Position',[2,1.2,0])
                case 3
                    lg = legend(flip(bPlot),{'Inertial';'Crossover';'Viscoelastic'},'Location','eastoutside','FontSize',12,'Position',[0.89,0.45,0.1,0.14]);
            end            
            text(2.6,.9,1,jointNames{jj},'Color','k','FontSize',14);
    end
    
    text(3.64,0.445,{'Dominant';'Forces'},'HorizontalAlignment','center')
    for jj = 1:3
        switch jj
            case 1
                set(aa(jj),'Position',[.0716,.1079,.1973,.815])
            case 2
                set(aa(jj),'Position',[.3759,.1079,.1973,.815])
            case 3
                set(aa(jj),'Position',[.6805,.1079,.1973,.815])
        end
    end
    
    miniData = mini_legend_info();
    if showStance
        miniLegSeq = [3,2,1,6,5,4,9,8,7];
    else
        miniLegSeq = 1:9;
    end
    for ii = 1:9
        ind = miniLegSeq(ii);   
        axes('pos',miniData{ii,1}); 
        [YourImage, ~, ImageAlpha] = imread(miniData{ind,2});
        h = imshow(YourImage);
        h.AlphaData = ImageAlpha;
    end
%% 1x3 - Data Only - Line Graphs Figure
modell2p = find(~all(isnan([phase_data(:,:,1),phase_data(:,:,2),phase_data(:,:,3)]),2));
cm = jet(length(modell2p)); jointNames = {'Hip';'Knee';'Ankle'};
figure('Position',[602,508,1318,481],'Name',strideName); subSize = []; bPlot = [];
    tiledlayout(1,3,'TileSpacing','compact')
    for jj = 1:3
        aa(jj) = subplot(1,3,jj);
            yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on
            bPlot = plot(phase_data(modell2p,:,jj)','-s','LineWidth',3);
            % Set scale colors
            for ii = 1:length(modell2p)
               bPlot(ii).Color = cm(ii,:); 
               bPlot(ii).MarkerFaceColor = cm(ii,:);
            end
            % Tweak tick info on axes
            set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]},'FontSize',10)
            % Add labels and title
            ylim([-1, 1]); grid on; pbaspect([1,1,1])
            set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
            set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
        switch jj
            case 1
                ylabel('Normalized Joint Torque','FontSize',16);
                set(aa(jj),'Position',[.0716,.1079,.1973,.815])
            case 2
                xlabel('Stride Phase','FontSize',16,'Position',[2,-1.3,-1]);
                set(aa(jj),'Position',[.3433,.1079,.1973,.815])
                title(['Data Joint Torques - ',strideName],'FontSize',16,'Position',[2,1.2,0])
            case 3
                set(aa(jj),'Position',[.6168,.1079,.1973,.815])
                lg = legend(flip(bPlot),flip(dTorque_legender(modell2p)),'Location','eastoutside','Orientation','vertical','FontSize',12);
        end
        text(2.6,.9,1,jointNames{jj},'Color','k','FontSize',14);
        
    end
    
    for jj = 1:3
        switch jj
            case 1
                set(aa(jj),'Position',[.0716,.1079,.1973,.815])
            case 2
                set(aa(jj),'Position',[.3433,.1079,.1973,.815])
            case 3
                set(aa(jj),'Position',[.6168,.1079,.1973,.815])
        end
        set(lg,'Position',[0.8408,.2786,.151,.5])
    end
    %legend(flip(bPlot),flip(cellfun(@num2str,num2cell(round(lengthVals(modell2p),2)),'UniformOutput',false)),'Location','eastoutside');
%% HIP ONLY, Cross Curves and Simple Lines
    range = 8; joint = 1; jointNames = {'Hip';'Knee';'Ankle'}; cm = [0 0 1; .5 0 .5; 1 0 0];
    titleSize= 20; labelSize = 16; inWordSize = 16; lineWidth = 3; tickSize = 14;
    if ~showStance
        strideName = 'Swing';
        range2plot = 1:37;
        crossOffset = .2425;
    else
        strideName = 'Stance';
        range2plot = 38:100;
        crossOffset = 1.03;
    end
    figure('Position',[158,313,1761,666]);
    aa(1) = subplot(1,2,1);
        for len = 1:numLens
            %r2p = range_finder(range,len,rangeMat);
            r2p = rangeMat(len,1):rangeMat(len,3);
            corrIn(joint,len) = corr(torques_active{len}(r2p,joint),torques_inertial{len}(r2p,joint));
            corrMs(joint,len) = corr(torques_active{len}(r2p,joint),torques_muscle{len}(r2p,joint));
            corrGv(joint,len) = corr(torques_active{len}(r2p,joint),torques_grav{len}(r2p,joint));
            if range >=5
                corrLd(joint,len) = corr(torques_active{len}(r2p,joint),torques_load{len}(r2p,joint));
            end
        end
        [x0,y0] = find_crossover_point(corrMs(1,:),corrIn(1,:),lengthVals);
        semilogx([x0 x0],[y0-.5 y0+.5],'Color',[.75,.75,.75],'LineWidth',lineWidth); hold on; grid on;
        semilogx(lengthVals,corrIn(joint,:),'r','LineWidth',lineWidth);
        semilogx(lengthVals,corrMs(joint,:),'b','LineWidth',lineWidth);
        ylim([-1.2,1.2]); ylabel('Relative Contribution','FontSize',labelSize); xlabel('x Size Relative to Rat','FontSize',labelSize);
        set(gca,'YTick',-1:.5:1,'FontSize',tickSize)
        xlim([lengthVals(1),lengthVals(end)]);
        text(0.06089,.8725,1,'Viscoelastic Forces','Color','b','FontSize',inWordSize);
%         text(0.0727,-0.8526,1,'Inertial Forces','Color','r','FontSize',inWordSize);
        text(0.06089,-0.702,1,'Inertial Forces','Color','r','FontSize',inWordSize);
        text(0.03,1.43,1,'A','FontSize',25)
    % Add animal labels
        xline(.375,'g-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
        xline(1,'c-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
        xline(2.63,'m-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
        xline(12,'y-.','LineWidth',2,'HandleVisibility','off','Alpha',.6)
        text(.1316,-1.1123,1,'Stick Insect \rightarrow','FontSize',15); 
        text(.5722,-1.1123,1,'Rat \rightarrow','FontSize',inWordSize); 
        text(1.4349,-1.1123,1,'Cat \rightarrow','FontSize',inWordSize); 
        text(5.2944,-1.1123,1,'Horse \rightarrow','FontSize',inWordSize); 
    % Add crossover info
        text(1.2*x0,y0,num2str(round(x0,2)),'FontSize',inWordSize)
        title(['Relative Contributions of Consituent Torques - Hip ',strideName],'FontSize',titleSize,'Units','normalized','Position',[.5,1.06,0])
        set(gca,'Position',[0.049,0.23,0.44,0.68])
    aa(2) = subplot(1,2,2,'FontSize',tickSize);
        [~,loc1] = min(abs(lengthVals-.05)); [~,loc2] = min(abs(lengthVals-x0)); [~,loc3] = min(abs(lengthVals-12));  
        modell2p = [loc1,loc2,loc3];
        yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on
        bPlot = plot(linspace(0,100,length(range2plot)),torques_act{1}(range2plot,modell2p),'LineWidth',lineWidth);
        %bPlot = plot(phase_model(modell2p,:,1)','-s','LineWidth',3);
        % Set scale colors
        for ii = 1:length(modell2p)
           bPlot(ii).Color = cm(ii,:); 
           bPlot(ii).MarkerFaceColor = cm(ii,:);
        end
        % Tweak tick info on axes
        %set(gca,'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]})
        %set(gca,'XTick',[1,2,3],'XTickLabel',{['Early ',strideName];['Mid ',strideName];['Late ',strideName]},'FontSize',12)
        % Add labels and title
        ylim([-1, 1]); grid on; pbaspect([1,1,1])
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
        ylabel('Normalized Joint Torque','FontSize',labelSize);
        xlabel([strideName,' (%)'],'FontSize',labelSize,'Position',[50,-1.5,-1]);
        title(['Modeled Active Muscle Torque - Hip ',strideName],'FontSize',titleSize,'Units','normalized','Position',[.5,1.06,0])
        lg = legend(flip(bPlot),{'Inertial';'Crossover';'Viscoelastic'},'Location','eastoutside','FontSize',inWordSize,'Position',[0.8703,0.5579,0.1,0.1]);
        text(127,0.42,{'Dominant';'Forces'},'HorizontalAlignment','center','FontSize',labelSize)
        text(-22.6,1.2,1,'B','FontSize',25)
        set(gca,'Position',[0.56,0.23,0.32,0.68])
        
        miniData = mini_legend_info_hip();
        if range == 8
            miniLegSeq = [3,2,1,6,5,4,9,8,7];
        else
            miniLegSeq = 1:9;
        end
        for ii = 1:3
            ind = miniLegSeq(ii);   
            axes('pos',miniData{ii,1}); 
            [YourImage, ~, ImageAlpha] = imread(miniData{ind,2});
            h = imshow(YourImage);
            h.AlphaData = ImageAlpha;
        end
%% HIP ONLY, Data compare
    joint = 1; jointNames = {'Hip';'Knee';'Ankle'}; modell2p = 1:numLens; cm = turbo(length(lengthVals));
    titleSize= 20; labelSize = 16; inWordSize = 16; lineWidth = 3; tickSize = 14;
    subPlotPos = {[0.0716,0.58,0.34,0.34];[0.4927,0.58,0.34,0.34];...
                  [0.0716,0.11,0.34,0.34];[0.4927,0.11,0.34,0.34]};
    legPos = {[];[0.8419,0.7,0.1247,0.09];[];[0.8419,0.18,0.155,0.198]};
    lineStyles = {'-';':';'-';'--';'--';'-.';':';'-';'-';'-';'-'};
    strideNames = {'Swing';'Swing';'Stance';'Stance'};
    range2plot = {1:37; 1:37; 38:100; 38:100};
    figure('Position',[446,2,1475,994]);
    for jj = 1:4
        if mod(jj,2) == 0 
            aa(jj) = subplot(2,2,jj,'FontSize',tickSize,'Position',subPlotPos{jj});
                d = squeeze(torqueNorm(range2plot{jj},joint,:));
                data2plot = ~all(isnan(d)); data2plot_inds = find(data2plot);
                aPlot = plot(linspace(0,100,length(d)),d(:,data2plot),'LineWidth',lineWidth); grid on; pbaspect([1.5,1,1])
                    % Set scale colors
                    for ii = 1:length(find(data2plot))
                       dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{data2plot_inds(ii)}),'',' ('))),2};
                       [~,minVal] = min(abs(dataScale-lengthVals));
                       dataInds(ii) = minVal;
                       aPlot(ii).Color = cm(minVal,:); 
                       aPlot(ii).MarkerFaceColor = cm(minVal,:); 
                       aPlot(ii).LineStyle = lineStyles{data2plot_inds(ii)};
                    end
                set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'},'FontSize',tickSize)
                title(['Active Muscle Torque from Data - ',strideNames{jj}],'FontSize',titleSize,'Units','normalized','Position',[.5,1.06,0]);
                legend(flip(aPlot),flip(dTorque_legender(data2plot_inds)),'FontSize',12,'Location','eastoutside',...
                    'Units','normalized','Position',legPos{jj});
                ylabel('Normalized Joint Torque','FontSize',labelSize,'Units','normalized','Position',[-0.08,0.5,0]);
                xlabel([strideNames{jj},' (%)'],'FontSize',labelSize);
        else
            aa(jj) = subplot(2,2,jj,'FontSize',tickSize,'Position',subPlotPos{jj});
                yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on;
                d = squeeze(torqueNorm(range2plot{jj+1},joint,:));
                data2plot = ~all(isnan(d)); data2plot_inds = find(data2plot);
                    for ii = 1:length(find(data2plot))
                       dataScale = scaleLookup{contains(scaleLookup(:,1),char(extractBetween(string(dTorque_legender{data2plot_inds(ii)}),'',' ('))),2};
                       [~,minVal] = min(abs(dataScale-lengthVals));
                       dataInds(ii) = minVal; 
                    end
                model2plot = unique(dataInds);
                bPlot = plot(linspace(0,100,length(range2plot{jj})),torques_act{joint}(range2plot{jj},model2plot),'LineWidth',lineWidth); grid on; pbaspect([1.5,1,1])
                % Set scale colors
                for ii = 1:length(model2plot)
                   bPlot(ii).Color = cm(model2plot(ii),:); 
                   bPlot(ii).MarkerFaceColor = cm(model2plot(ii),:);
                end
                % Add labels and title
                set(gca,'YLim',[-1 1],'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
                ylabel('Normalized Joint Torque','FontSize',labelSize,'Units','normalized','Position',[-0.08,0.5,0]);
                xlabel([strideNames{jj},' (%)'],'FontSize',labelSize);
                title(['Modeled Active Muscle Torque - ',jointNames{joint},' ',strideNames{jj}],'FontSize',titleSize,'Units','normalized','Position',[.5,1.06,0])
        end
    end
%% HIP ONLY, Predicted
modell2p = [26,13,1]; cm = hsv(length(modell2p));
titleSize= 20; labelSize = 16; inWordSize = 16; lineWidth = 3; tickSize = 14;
figure('Position',[171,33,1739,946]);
for ii = 1:length(modell2p)
    plot(torques_full{1}(:,modell2p(ii)),'Color',cm(ii,:),'LineWidth',lineWidth); hold on
end
grid on; pbaspect([3,1,1])
xline(37,'m--','LineWidth',2)
set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'},'FontSize',tickSize)
legend({'Brontosaurus';'Rat';'Fruit Fly'},'Location','eastoutside')
xlabel('Stride (%)','FontSize',labelSize)
ylabel('Normalized Joint Torque','FontSize',labelSize,'Units','normalized','Position',[-0.06,0.5,0]);
title('Predicted Active Joint Torque','FontSize',titleSize)
%% HIP ONLY, Motion Phase
modell2p = [1];
figure;
subplot(2,2,1)
    dInd = 8;   % stick insect data ind
    range2plot = rangeMat(dInd,1):rangeMat(dInd,3);
    jm = objCell{dInd}.theta_motion(range2plot,1); jm = interp1(1:length(jm),jm,linspace(1,length(jm),100))'; tm = torques_full{1}(:,dInd);
    plot(1:length(jm),[jm tm],'LineWidth',3)
    legend({'Joint Motion';'Torque'})
subplot(2,2,3)
    dInd = 22;   % horse data ind
    range2plot = rangeMat(dInd,1):rangeMat(dInd,3);
    jm = objCell{dInd}.theta_motion(range2plot,1); jm = interp1(1:length(jm),jm,linspace(1,length(jm),100))'; tm = torques_full{1}(:,dInd);
    plot(1:length(jm),[jm tm],'LineWidth',3)
    legend({'Joint Motion';'Torque'})
subplot(2,2,2)
    for ii = 1:numLens
        range2plot = rangeMat(ii,1):rangeMat(ii,3);
        jm = objCell{ii}.theta_motion(range2plot,1); jm = interp1(1:length(jm),jm,linspace(1,length(jm),100))'; tm = torques_full{1}(:,ii);
        %%phaseDiffs(ii) = phdiffmeasure(jm, tm)*(180/pi);
        phaseDiffs(ii) = corr(jm,tm);
    end
    semilogx(lengthVals, phaseDiffs,'LineWidth',3)
    legend({'Phase Diffs'})
%% Compare Model Animals at the right SIZE and WALK profile to Data Animals
jointNames = {'Hip';'Knee';'Ankle'};
cm = [1,0,0;.5,0,0;0,1,0;0,.5,0;0,0,1;0,0,.5];
figure('Position',[602,508,1318,481],'Name','Stance'); subSize = []; bPlot = [];
    %tiledlayout(3,1,'TileSpacing','compact')
    stanceMeaner = @(inMat) [mean(inMat(38:58,:));mean(inMat(59:80,:));mean(inMat(81:end,:))];
    rat = [stanceMeaner(torqueNorm(:,:,2));stanceMeaner(objTorques_tot{1})];
    cat = [stanceMeaner(torqueNorm(:,:,5));stanceMeaner(objTorques_tot{5})];
    horse = [stanceMeaner(torqueNorm(:,:,11));stanceMeaner(objTorques_tot{9})];
    for jj = 1:3
        aa(jj) = subplot(1,3,jj);
            %modell2p = 1:numLens; 
            yline(0,'LineWidth',2,'Alpha',.1,'HandleVisibility','off'); hold on
            bPlot = plot([rat(1:3,jj),rat(4:6,jj),cat(1:3,jj),cat(4:6,jj),horse(1:3,jj),horse(4:6,jj)],'-s','LineWidth',3);
%             % Set scale colors
            for ii = 1:6
               bPlot(ii).Color = cm(ii,:); 
               bPlot(ii).MarkerFaceColor = cm(ii,:);
            end
            % Tweak tick info on axes
            set(gca,'XTick',[1,2,3],'XTickLabel',{['Early Stance'];['Mid Stance'];['Late Stance']})
            %set(gca,'XTickLabel',{['Early Stance'];['Mid Stance'];['Late Stance']})
            % Add labels and title
            ylim([-1, 1]); grid on; pbaspect([1,1,1])
            set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
            text(1.03,.9,1,jointNames{jj},'Color','k','FontSize',14);
        if jj == 1
            ylabel('Normalized Joint Torque','FontSize',14);
        end
        if jj == 2
            xlabel('Stride Phase','FontSize',14,'Position',[2,-1.3,-1])
        end
        text(1.03,.9,1,jointNames{jj},'Color','k','FontSize',14);
        if jj == 2
            title(['Active Joint Torques - Stance'],'FontSize',16,'Position',[2,1.2,0])
        end
    end
    lg = legend(aa(2), {'Rat Data';'Rat Model';'Cat Data';'Cat Model';'Horse Data';'Horse Model'},'Location','eastoutside',...
        'Position',[.3,.028,.4347,.0478],'Orientation','horizontal','FontSize',14);
%% Compare Model to Data, focus on one joint
temp = [1,5,9]; temp2 = [2,5,11]; animals = {'Rat';'Cat';'Horse'};
figure('Position',[762,-7,1159,1003])
for ii = 1:3
    subplot(3,1,ii)
    plot(linspace(0,100,63),objTorques_tot{temp(ii)}(38:end,1),'LineWidth',3);hold on;plot(linspace(0,100,63),torqueNorm(38:end,1,temp2(ii)),'LineWidth',3)
    legend({'Model';'Data'});ylabel('Normalized Joint Torque');
        xlabel('Stance %')
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'})
    title([animals{ii},' walking like ',animals{ii}])
end
%% Compare Active Torques Across Walking Types WAVEFORMS
figure('Position',[1,1,1920,1003],'Name','test'); c2p = {[1,2,3];[4,5,6];[7,8,9]}; jointNames = {'Hip';'Knee';'Ankle'}; animNames = {'Rat';'Cat';'Horse'};
counter = 1; waveColors = {'r';'g';'b'};
for animal = 1:3
    %animal = 3; % 1 = Rat, 2 = Cat, 3 = Horse
    hip = [objTorques_act{c2p{animal}(1)}(:,1),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
    kne = [objTorques_act{c2p{animal}(1)}(:,2),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
    ank = [objTorques_act{c2p{animal}(1)}(:,3),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
    
    for joint = 1:3
        aa = subplot(3,3,counter);
        
        jDat = [objTorques_act{c2p{animal}(1)}(:,joint),objTorques_act{c2p{animal}(2)}(:,joint),objTorques_act{c2p{animal}(3)}(:,joint)];
        plot(jDat,'LineWidth',3); grid on
        colororder(aa,[1 0 0; 0 1 0; 0 0 1]);
        title(jointNames{joint},'FontSize',20)
        set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'},'FontSize',12)
        xlabel('Stride %','FontSize',14); yline(0,'HandleVisibility','off');xline(37,'--m','HandleVisibility','off')
        if joint == 1
            ylabel({[animNames{animal},' Size'];'Normalized Joint Torque'},'FontSize',16,'FontWeight','bold')
        else
            ylabel('Normalized Joint Torque','FontSize',14)
        end
        if joint == 3
            legend({'Rat Walking';'Cat Walking';'Horse Walking'},'FontSize',12,'Location','eastoutside')
        end
        temp = get(aa,'Position');
        set(aa,'Position',[temp(1:2),.2,.2])
        counter = counter + 1;
    end
end
%% Compare Active Torques Across Walking Types LINES
animal = 1; % 1 = Rat, 2 = Cat, 3 = Horse
c2p = {[1,2,3];[4,5,6];[7,8,9]}; jointNames = {'Hip';'Knee';'Ankle'}; animNames = {'Rat';'Cat';'Horse'};
hip = [objTorques_act{c2p{animal}(1)}(:,1),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
kne = [objTorques_act{c2p{animal}(1)}(:,2),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
ank = [objTorques_act{c2p{animal}(1)}(:,3),objTorques_act{c2p{animal}(2)}(:,1),objTorques_act{c2p{animal}(3)}(:,1)];
figure('Position',[762,-7,1159,1003],'Name','test'); counter = 1;
for joint = 1:3
    subplot(3,1,joint)
    jDat = [objTorques_act{c2p{animal}(1)}(:,joint),objTorques_act{c2p{animal}(2)}(:,joint),objTorques_act{c2p{animal}(3)}(:,joint)];
    jDat = [mean(jDat(1:12,:));mean(jDat(13:24,:));mean(jDat(25:36,:))];
    plot(jDat,'-s','LineWidth',3); grid on
    legend({'Rat Walking';'Cat Walking';'Horse Walking'},'FontSize',12,'Location','eastoutside')
    title([jointNames{joint},' Active Torques - ',animNames{animal},' Size'],'FontSize',20)
    set(gca,'YTick',-1:.5:1,'YTickLabel',{'FLX  -1';'-.5';'0';'.5';'EXT  1'},'FontSize',12)
    set(gca,'XTick',[1,2,3],'XTickLabel',{['Early Swing'];['Mid Swing'];['Late Swing']})
    ylabel('Normalized Joint Torque','FontSize',14)
end
%% Plot Poly Waveforms for Model
figure('Position',[962,2,958,994],'Name',strideName); jointNames = {'Hip';'Knee';'Ankle'}; cm = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];
curvars = whos;
if any(contains({curvars.name},'crossovers'))
    if showStance == 1
        cross = crossovers(8,:);
    else
        cross = crossovers(4,:);
    end
else
    if showStance == 1
        cross = [.62,.22,.54];
    else
        cross = [.75,.36,1.3];
    end
end
outSlopesModel = find_all_slopes_poly(phase_model);
[outSlopesData,gofData] = find_all_slopes_poly(phase_data);
for ii = 1:3
    subplot(3,1,ii)
    semilogx(lengthVals,outSlopesModel(:,:,ii),'LineWidth',2);hold on;
    yline(0); ylim([-1 1])
    xline(cross(ii))
    % Overlay data
    for jj = 1:size(outSlopesData,1)
        dScale = dTorque_scaler(jj);
        if ~all(outSlopesData(jj,:,1)==0)
            scatter([dScale,dScale,dScale],outSlopesData(jj,:,ii),35,cm,'filled')
        end
    end
    title([jointNames{ii},' - ',strideName],'Fontsize',16);
    legend({'x^2','x','c'},'Location','northwest','FontSize',10)
end
%% Function: Mini Legend Locations
function outData = mini_legend_info_hip()
    rootDir = 'G:\My Drive\Rat\Swing Paper\MiniLegend\'; axW = .1; axH = .1; y = .08;
    outData = {[.535, y, axW, axH],[rootDir,'Hip-1.png'];...
               [.675, y, axW, axH],[rootDir,'Hip-2.png'];...
               [.8, y, axW, axH],[rootDir,'Hip-3.png']};
end
%% Function: Mini Legend Locations
function outData = mini_legend_info()
    rootDir = 'G:\My Drive\Rat\Swing Paper\MiniLegend\'; axW = .1; axH = .1; y = .08;
    outData = {[.02, y, axW, axH],[rootDir,'Hip-1.png'];...
               [.12, y, axW, axH],[rootDir,'Hip-2.png'];...
               [.22, y, axW, axH],[rootDir,'Hip-3.png'];...
               [.322, y, axW, axH],[rootDir,'Knee-1.png'];...
               [.425, y, axW, axH],[rootDir,'Knee-2.png'];...
               [.525, y, axW, axH],[rootDir,'Knee-3.png'];...
               [.625, y, axW, axH],[rootDir,'Ankle-1.png'];...
               [.728, y, axW, axH],[rootDir,'Ankle-2.png'];...
               [.828, y, axW, axH],[rootDir,'Ankle-3.png']};
end
%% Find slope POLY
function [outCoeffs,gof] = find_line_slope_poly(inLine)
    if any(isnan(inLine))
        outCoeffs = [0,0,0]; gof.rsquare = -2;
        return
    end
    [xData, yData] = prepareCurveData( [], inLine );
    % Set up fittype and options.
    ft = fittype( 'poly2' );
    % Fit model to data.
    [fitresult,gof] = fit( xData, yData, ft ,'Normalize', 'on' );
    outCoeffs = coeffvalues(fitresult);
end
%% Find slopes of whole data chunk POLY
function [outSlopes,gofOut] = find_all_slopes_poly(inPhase)
    for ii = 1:size(inPhase,1)
        for jj = 1:3
            [outSlopes(ii,:,jj),gof] = find_line_slope_poly(inPhase(ii,:,jj));
            gofOut(ii,jj) = gof.rsquare;
        end
    end
end
%% Range_finder: Output Indices for Range of Motion
function range2plot = range_finder(range,len,rangeMat)
    beg = rangeMat(len,1); mid = rangeMat(len,2); ennd = rangeMat(len,3);
    switch range
        case 1
            endInd = beg + floor((mid-beg)/3);
            range2plot = beg:endInd;
        case 2
            endInd1 = beg + floor((mid-beg)/3);
            endInd2 = beg + floor(2*(mid-beg)/3);
            range2plot = endInd1:endInd2;
        case 3
            endInd = beg + floor(2*(mid-beg)/3);
            range2plot = endInd:mid;
        case 4
            range2plot = rangeMat(len,1):rangeMat(len,2);
        case 5
            endInd = mid + floor((ennd-mid)/3);
            range2plot = mid:endInd;
        case 6
        case 7

        case 8
            range2plot = rangeMat(len,2):rangeMat(len,3);
    end
end
%% Find_Crossover_Point
function [outX,outY] = find_crossover_point(corr1,corr2,lengthVals)
    % Scan across the two correlation vectors and find the point at which they intersect
    aa = corr1; bb = corr2;    
    for ii = 2:length(aa)
        if aa(ii) < bb(ii) && aa(ii-1) > bb(ii-1) % if we just passed the crossover point, find the intersection
            m2 = (bb(ii)-bb(ii-1))/(lengthVals(ii)-lengthVals(ii-1));
            x = lengthVals(ii-1);
            y2 = bb(ii-1);
            m1 = (aa(ii)-aa(ii-1))/(lengthVals(ii)-lengthVals(ii-1));
            y1 = aa(ii-1);
            outX = ((m2-m1)*x+y1-y2)/(m2-m1);
            outY = (m2*y1-m1*y2)/(m2-m1);
        end
    end
end