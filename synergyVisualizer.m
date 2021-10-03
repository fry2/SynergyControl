function [colorlog,widthlog] = synergyVisualizer(obj,W,H,to_save)
    % All forces must be positive for this to work. Line widths are based on activation level, can't be negative.\
    scrSz = get(groot, 'ScreenSize');
    scW = scrSz(3);
    [beg,ennd,~] = find_step_indices(obj);
    div = 500;
    xx = floor(linspace(beg,ennd,div));
    tempcolor = hot(256);
    
            limrecorder =zeros(length(xx),6);
            colorlog = [];
            widthlog = [];
    %                 v = VideoWriter(['G:\My Drive\Rat\ForceMinimizationAnalysis\Videos\',mintypes{minmethod},'_',datestr(datetime('now'),'yyyymmdd')],'Motion JPEG AVI');
    %                 v.Quality = 100;
    %                 v.FrameRate = 5;
    %                 open(v);
    %% Check all open figures to see if planned figures are already open. Close all currently existing planned figures.
    figHandles = get(groot, 'Children');
    if ~isempty(figHandles)
        priorFigures = contains({figHandles(:).Name},{'ForceFig','SynergyFig'});
        close(figHandles(priorFigures))
    end

    forces = W*H;
    
    forceOrSyn = 2;
    switch forceOrSyn
        case 1
            displayFig = figure('color','white','Position',[-(scW-10)/2 130 scW/2.03 985],'name','ForceFig');
            trim = floor(.2*length(tempcolor));
            tempcolor = tempcolor(trim:(end-trim),:);
            CM = colormap(flipud(tempcolor));
            colorMat = CM;
            maxforce = max(forces,[],'all');
        case 2
            % choose this to display synergy coloring
            displayFig = figure('color','white','Position',[-(scW-10) 130 scW/2.03 985],'name','SynergyFig');
            relW = W./max(W);
            bigH = (H'.*max(W))';
            relH = (bigH./max(bigH')');
            cm = colormap(jet(3*size(W,2)));
            cm = cm(3:3:end,:);
    end
    
    for timecount = 1:2:length(xx)
        if forceOrSyn == 2
            muscleweights = relW.*relH(:,timecount)';
            muscleweights = muscleweights./max(muscleweights,[],2);
            synergycolors = zeros(38,3);
            for jj = 1:38
                holder = muscleweights(jj,:)'.*cm;
                if any(mean(holder)==0)
                    emptyCol = find(mean(holder)==0);
                    [ii,~,v] = find(holder');
                    temp = accumarray(ii,v,[],@mean)';
                    for ind = 1:length(emptyCol)
                        temp(emptyCol(ind)) = 0;
                    end
                    synergycolors(jj,:) = temp;
                else
                    [ii,~,v] = find(holder');
                    synergycolors(jj,:) = accumarray(ii,v,[],@mean)';
                end
            end
            colorMat = synergycolors;
        else
            colorvec = 1+floor(length(CM)*(forces(:,timecount)./maxforce));
            colorvec(colorvec(:,1)==156,:)=155;
            %colorvec = floor(rescale((forces(:,timecount)),1,length(CM)));
            %widthvec = rescale(colorvec,.1,3);
            widthvec = 2*(forces(:,timecount));
            colorlog(:,timecount) = colorvec;
            widthlog(:,timecount) = widthvec;
        end
        widthvec = 2*(forces(:,timecount));
        figure(displayFig)
        cla(displayFig.CurrentAxes)
        %% Plot muscle paths
        for musnum = 1:38
            muscle = obj.musc_obj{musnum};
            musclemat =zeros(size(muscle.pos_attachments,1),3);
            for hh = 1:size(muscle.pos_attachments,1)
                temp = 1000*muscle.pos_attachments{hh,4}(xx(timecount),:);
                musclemat(hh,:) = temp;
            end

            if forceOrSyn == 1
                colorSelector = colorvec(musnum,1);
            else
                colorSelector = musnum;
            end

            plot3(musclemat(:,1),musclemat(:,2),musclemat(:,3),'Color',colorMat(colorSelector,:),'Linewidth',widthvec(musnum,1))
            hold on
        end
        %% Plot joints and joint connections
        jointmat = zeros(4,3);
        for jointnum = 1:3
            joint = obj.joint_obj{jointnum};
            temp  = 1000*joint.sim_position_profile(xx(timecount),:);
            jointmat(jointnum,:) = temp;
        end
        jointmat(4,:) = 1000*obj.musc_obj{20}.pos_attachments{end,4}(xx(timecount),:);
        plot3(jointmat(:,1),jointmat(:,2),jointmat(:,3),'k--','LineWidth',2,'Marker','o','MarkerFaceColor','k','MarkerSize',5)
        hold on
        %% Prepare axes proportions for proper leg dimensions
            limrecorder(timecount,:) = [xlim,ylim,zlim];
            xlimits = [-60 60];
            ylimits = [-80 20];
            zlimits = [0 30];
            xlim(xlimits);
            ylim(ylimits);
            zlim(zlimits);
            a =sum(abs(xlimits));
            b =sum(abs(ylimits));
            c =sum(abs(zlimits));
            axeslengths = [a b c];
            normedaxes = axeslengths/norm(axeslengths);
            pbaspect(normedaxes)
            view([0 90])
        title('Synergy Activation During Walking','FontSize',15)
        %cbar = colorbar('Ticks',[0;1],'TickLabels',{'low act','high act'},'FontSize',15);
        grid on
        pause(.001)
        if to_save
            frame = getframe(displayFig);
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if timecount == 1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',.1); 
            else 
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',.1); 
            end 
        end
%                 writeVideo(v,frame);
        hold off
    end
%             close(v)
    clear forces v
end