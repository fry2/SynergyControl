function plotWH(inArray,W,H,saver)
    % Plot the synergy components W and H
    
    scrSz = get(groot, 'ScreenSize');
    scW = scrSz(3);
    
    [r,c] = size(inArray);
    
    if r>c
        inArray = inArray';
    end
    
    if iscell(W)
        Wth = W{2}*ones(c,1);
        W = W{1};
    else
        Wth = 0;
    end
    
    
    isforce = abs(max(W*H,[],'all')) > 1;
    
    if isforce
        relW = W./max(W);
        bigH = (H'.*max(W))';
    else
        relW = W;
        bigH = H;
    end
    
    recompiled = (W*H)';
    
    figHandles = get(groot, 'Children');
    if ~isempty(figHandles)
        priorFigures = contains({figHandles(:).Name},{'SynFig','CFig','RecombFig'});
        close(figHandles(priorFigures))
    end
    axisticker = 1:c;
    synfig = figure('color','white','Position',[-(scW-10) 130 scW/2 985],'name','SynFig');
    cm = colormap(jet(3*size(W,2)));
    for i = 1:size(W,2)
        holder = relW(:,i);
        subplot(size(relW,2),1,i)
        bar(holder,'FaceColor',cm(i*3,:))
        if size(Wth,1) ~= 1
            hold on
            plot(axisticker,Wth,'r--','LineWidth',1.2)
        end
        set(gca,'XTick',axisticker,'TickLength',[.003 .003])
        if i == 1
            title('Relative Activation of Individual Muscles','FontSize',16)
            if size(Wth,1) ~= 1
                title(['Relative Activation of Individual Muscles. Threshold: ',num2str(round(Wth(1),2))],'FontSize',16)
            else
                title('Relative Activation of Individual Muscles','FontSize',16)
            end
        elseif i == size(W,2)
            xlabel('Muscle #','FontSize',14)
        end
        xlim([0 size(W,1)+1])
        ylabel(num2str(i))
    end
    
    cfig = figure('color','white','Position',[-(scW-10)/2 130 scW/2 985],'name','CFig');
    yLims = [min(bigH,[],'all') max(bigH,[],'all')];
    for j = 1:size(bigH,1)
        holder = bigH(j,:);
        subplot(size(bigH,1),1,j)
        plot(linspace(0,100,size(bigH,2)),holder,'LineWidth',3,'Color',cm(j*3,:))
        ylim(yLims)
        ylabel(num2str(j))
        if j == 1
            title('Synergy Activation During Stride','FontSize',16)
            if size(Wth,1) ~= 1
                title(['Synergy Activation During Stride. Threshold: ',num2str(round(Wth(1),2))],'FontSize',16)
            else
                title('Synergy Activation During Stride','FontSize',16)
            end
        elseif j == size(bigH,1)
            xlabel('% Stride','FontSize',14)
        end
    end


    recombfig = figure('color','white','Position',[50 50 800 700],'name','RecombFig');
    yMax = max([max(inArray,[],'all') max(recompiled,[],'all')]);
    yMin = min([min(inArray,[],'all') min(recompiled,[],'all')]);
    subplot(3,1,1)
        plot(linspace(0,100,size(inArray,2)),inArray')
        ylim([yMin yMax])
        title('Original Signal')
    subplot(3,1,2)
        plot(linspace(0,100,size(inArray,2)),recompiled)
        ylim([yMin yMax])
        title('Recompiled Signal from NNMF Components')
    subplot(3,1,3)
        differ = mean(abs((inArray'-recompiled)));
        bar(differ)
%         differ = diag(corr(recompiled,inArray'));
%         bar(differ)
        %plot(linspace(0,100,size(inArray,2)),differ)
        if size(Wth,1) ~= 1
            title(['Mean Difference Between Signals. Threshold: ',num2str(round(Wth(1),2))])
        else
            title('Mean Difference Between Signals')
        end
    
    if saver
        saveas(cfig,['G:\My Drive\Rat\SynergyControl\OutputFigures\Images\',datestr(datetime('now'),'yyyymmdd'),'_','actfig.png']);
        saveas(synfig,['G:\My Drive\Rat\SynergyControl\OutputFigures\Images\',datestr(datetime('now'),'yyyymmdd'),'_','wgtfig.png']);
        saveas(recombfig,['G:\My Drive\Rat\SynergyControl\OutputFigures\Images\',datestr(datetime('now'),'yyyymmdd'),'_','recfig.png']);
    end
end