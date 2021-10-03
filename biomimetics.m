function biomimetics(obj)
    % BFA 33
    % Pect 5
    % Semimembranosus 31
    % VastusIntermedius 18
    % ExtDiGLong 20
    % MedGastroc 17
    % TibialisAnt 22
    sensemuscname = 'BFA';
    sensecasenum = 1;
    joint = 1;
    axis = 1;
    [beg,ennd,~] = find_step_indices(obj);
    %Higher div, the longer the program takes.
    div = 100;
    xx = floor(linspace(beg,ennd,div));
    xx = obj.sampling_vector;
    if 0
        muscles = [11 11 33 33 27 27;...
             1 2 1 2 2 3;...
             2 1 2 1 3 2];
         limitsy = [9.5 10 3 .8 4.6 0;...
             16 16.5 8 4.2 5.6 6];
        musclenames = {'Biceps Femoris Posterior (Hip)';'Biceps Femoris Posterior (Knee)';'Biceps Femoris Anterior (Hip)';...
            'Biceps Femoris Anterior (Knee)';'Plantaris (Knee)';'Plantaris (Ankle)'};
        jointnames = {'Hip';'Knee';'Ankle'};
         for i=1:length(muscles)
            subplot(3,2,i)
            m = obj.musc_obj{muscles(1,i)};
            [moment_arm_profile,~] = compute_muscle_moment_arm(obj,m,axis,muscles(2,i),0,0);
            momprof = moment_arm_profile(:,3);
            %momprof = (momprof-min(momprof))/(max(momprof)-min(momprof));
            bb = m.muscle_length_profile(xx);
            %bb = (bb-min(bb))/(max(bb)-min(bb));
            tdlo = [];
            tdlo(1,1) = 1000*bb(6);
            tdlo(1,2) = momprof(6);
            tdlo(2,1) = 1000*bb(51);
            tdlo(2,2) = momprof(51);
            plot(1000*bb,momprof,'k','LineWidth',1)
            hold on
            %Touch Down Marker
            td = plot(tdlo(1,1),tdlo(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0]);
            %Lift Off Marker
            lo = plot(tdlo(2,1),tdlo(2,2),'go','MarkerSize',10,'MarkerFaceColor',[0 1 0]);
            cc = [1000*bb,momprof];
            ee = cc(10:20:end,:);
            ff = cc(11:20:end,:);
            C = ff([1;1]*(1:size(ff,1)),:);
            C(1:2:end,:) = ee;
                                cc2 = moment_arm_profile(:,2);
            ee2 = cc2(10:20:end,:);
            ff2 = cc2(11:20:end,:);
            C2 = ff2([1;1]*(1:size(ff2,1)),:);
            C2(1:2:end,:) = ee2;
            for j = 1:2:length(C)-1
%                         arrow(C(j,:),C(j+1,:),'Length',5,'Width',0,'TipAngle',30,'Color','b')
                arrow(C(j,:),C(j+1,:),'Length',5,'Width',.001,'Color','b')
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    arrow FIXLIMITS
                    lastwarn('');
                end
            end
            grid on
            title(musclenames{i})
            ylabel('Moment Arm Length (mm)')
            xlabel('Muscle Length (mm)')
            ylim(limitsy(:,i))
            switch i
                case 1
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','northeast')
                case 2
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','southeast')          
                case 3
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','southwest')
                case 4
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','southeast')
                case 5
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','northeast')
                case 6
                    legend([td,lo],{'Touch Down','Lift Off'},'Location','southeast')
            end
         end
         set(gcf,'Position',[800 10 1000 800])
         keyboard
    end
    %This generates the charles_compare_figure subplots, the main isolated joint figures. There's no sensitivity analysis in this one.
    if 0
        muscles = [33 5 31 18 17 22 8 9 36;...
                1 1 2 2 3 3 1 1 1];
        color = [1 0 0;0 1 0;0 0 1];
        limitsy = [-.4 -.2 -.3 0 -.15 -.05 -.143 -.063 -.08;...
                    .2 .2 0 .15 0 .15 .57 .043 .086];
        limitsx = [-30 -30 -145 -145 -50 -50 -30 -30 -30;...
            50 50 -45 -45 50 50 50 50 50];
        musclenames = {'Biceps Femoris Anterior';'Pectineus';'Semimembranosus';...
        'Vastus Intermedius';'Medial Gastrocnemius';'Tibialis Anterior';'Tensor Fascia Latae';'Obturator Externus';'Obturator Internus'};
        for i=1:length(muscles)
            if size(unique(obj.theta_motion(:,muscles(2,i))),1)==1
            else
                [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,i)},axis,muscles(2,i),0,0);
                figure
                if muscles(2,i) == 1
                    jointangle = -moment_arm_profile(:,2);
                    musclename = obj.joint_obj{muscles(2,i)}.name(4:end-1);
                    momentarm = -moment_arm_profile(:,3)/35.745;
                elseif muscles(2,i) == 2
                    jointangle = moment_arm_profile(:,2)-81.578;
                    musclename = obj.joint_obj{muscles(2,i)}.name(4:end);
                    momentarm = -moment_arm_profile(:,3)/35.745;
                else
                    jointangle = -moment_arm_profile(:,2)-20.0535;
                    musclename = obj.joint_obj{muscles(2,i)}.name(4:end-1);
                    momentarm = -moment_arm_profile(:,3)/35.745;
                end
                %momnorm = (momentarm-min(momentarm))/(max(momentarm)-min(momentarm));
                if muscles(2,i) == 2
                    plot(jointangle,-momentarm,'color',color(axis,:),'LineWidth',2)
                else
                    plot(jointangle,momentarm,'color',color(axis,:),'LineWidth',2)
                end
                grid off
                title(musclenames{i})
                ylabel('Moment Arm/thigh length')
                xlabel([musclename,' Angle (deg)'],'Interpreter','none')
                ylim(limitsy(:,i))
                xlim(limitsx(:,i))
                %saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper\complete_charles_comparisons\',obj.musc_obj{muscles(1,i)}.muscle_name(4:end),'_',datestr(datetime('now'),'yyyymmdd'),'.png'])
                %saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper\complete_charles_comparisons\',obj.musc_obj{muscles(1,i)}.muscle_name(4:end),'_',datestr(datetime('now'),'yyyymmdd'),'.eps'])
            end
        end
        keyboard
    end
    %Sensitivity analysis of the isolated joint moment arms
    if 0      
        muscles = [33 5 31 18 17 22 8 9 36;...
                1 1 2 2 3 3 1 1 1];
        color = [1 0 0;0 1 0;0 0 1];
        limitsy = [-.4 -.2 -.3 0 -.15 -.05 -.143 -.063 -.08;...
                    .2 .2 0 .15 0 .15 .57 .043 .086];
        limitsx = [-30 -30 -145 -145 -50 -50 -30 -30 -30;...
            50 50 -45 -45 50 50 50 50 50];
        musclenames = {'Biceps Femoris Anterior';'Pectineus';'Semimembranosus';...
        'Vastus Intermedius';'Medial Gastrocnemius';'Tibialis Anterior';'Tensor Fascia Latae';'Obturator Externus';'Obturator Internus'};
        switch sensemuscname
            case 'BFA'
                i = 1;
            case 'Pect'
                i = 2;
            case 'SM'
                i = 3;
            case 'VI'
                i = 4;
            case 'MG'
                i = 5;
            case 'TA'
                i = 6;
        end
        switch sensecasenum
            case 1
                linecolor = 'k';
                linetype = '-';
                pathpath = ['G:\My Drive\Rat\BiomimeticsPaper\Sensitivity\',sensemuscname,'\Base.png'];
            case 2
                linecolor = 'b';
                linetype = '-';
                pathpath = ['G:\My Drive\Rat\BiomimeticsPaper\Sensitivity\',sensemuscname,'\Origin+.png'];
            case 3
                linecolor = 'b';
                linetype = '--';
                pathpath = ['G:\My Drive\Rat\BiomimeticsPaper\Sensitivity\',sensemuscname,'\Origin-.png'];
            case 4
                linecolor = 'r';
                linetype = '-';
                pathpath = ['G:\My Drive\Rat\BiomimeticsPaper\Sensitivity\',sensemuscname,'\Ins+.png'];
            case 5
                linecolor = 'r';
                linetype = '--';
                pathpath = ['G:\My Drive\Rat\BiomimeticsPaper\Sensitivity\',sensemuscname,'\Ins-.png'];
        end
        if size(unique(obj.theta_motion(:,muscles(2,i))),1)==1
        else
            [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,i)},axis,muscles(2,i),0,0);
            figure
            if muscles(2,i) == 1
                jointangle = -moment_arm_profile(:,2);
                musclename = obj.joint_obj{muscles(2,i)}.name(4:end-1);
                momentarm = -moment_arm_profile(:,3)/35.745;
            elseif muscles(2,i) == 2
                jointangle = moment_arm_profile(:,2)-81.578;
                musclename = obj.joint_obj{muscles(2,i)}.name(4:end);
                momentarm = -moment_arm_profile(:,3)/35.745;
            else
                jointangle = -moment_arm_profile(:,2)-20.0535;
                musclename = obj.joint_obj{muscles(2,i)}.name(4:end-1);
                momentarm = -moment_arm_profile(:,3)/35.745;
            end
            %momnorm = (momentarm-min(momentarm))/(max(momentarm)-min(momentarm));
            if muscles(2,i) == 2
                plot(jointangle,-momentarm,'color',linecolor,'LineWidth',2,'LineStyle',linetype)
            else
                plot(jointangle,momentarm,'color',linecolor,'LineWidth',2,'LineStyle',linetype)
            end
            grid off
            title(musclenames{i})
            ylabel('Moment Arm/thigh length')
            xlabel([musclename,' Angle (deg)'],'Interpreter','none')
            ylim(limitsy(:,i))
            xlim(limitsx(:,i))
            saveas(gcf,pathpath)
        end
        return
    end
    %Generates the 6 simulation images (with the flyouts) for the BFP, BFA, and the Plantaris
    if 0
        muscles = [11 11 33 33 27 27;...
            1 2 1 2 2 3;...
            2 1 2 1 3 2];
        for ii = 1:size(muscles,2)
            compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,ii)},axis,muscles(2,ii),0,1);
        end
        keyboard
    end
    if 0
%                 muscles = [11 12 33 14 27 20;...
%                     1 1 1 1 2 2;...
%                     2 2 2 2 3 3];
%                 musclenames = {'Biceps Femoris Posterior';'Rectus Femoris';'Biceps Femoris Anterior';...
%                     'Semitendinosus Accessory';'Plantaris';'Extensor Digitorum Longus'};
        muscles = [11 11 33 33 27 27;...
             1 2 1 2 2 3;...
             2 1 2 1 3 2];
        zlimits = [9.8 9.8 1.4 1.4 .8 .8
            16 16 7.8 7.8 6 6];
        musclenames = {'Biceps Femoris Posterior (Hip)';'Biceps Femoris Posterior (Knee)';'Biceps Femoris Anterior (Hip)';...
            'Biceps Femoris Anterior (Knee)';'Plantaris (Knee)';'Plantaris (Ankle)'};
        jointnames = {'Hip';'Knee';'Ankle'};
        count2 = 1;
        for ii = 1:size(muscles,2)
            [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,ii)},axis,muscles(2,ii),0,0);
            jointangles = [-obj.joint_obj{1}.rec_angle_profile(xx)*(180/pi),...
                obj.joint_obj{2}.rec_angle_profile(xx)*(180/pi)-81.578,...
                -obj.joint_obj{3}.rec_angle_profile(xx)*(180/pi)-20.0535];
            %mnorm = (moment_arm_profile(:,3)-min(moment_arm_profile(:,3)))/(max(moment_arm_profile(:,3))-min(moment_arm_profile(:,3)));
            subplot(3,2,ii)
            zlim(zlimits(:,count2)')
            bb = fill3(jointangles(:,muscles(2,count2)),jointangles(:,muscles(3,count2)),moment_arm_profile(:,3),'r');
            %bb = fill3(jointangles(:,muscles(3,count2)),jointangles(:,muscles(2,count2)),mnorm,'r');
            hold on
            bb.FaceAlpha = .5;
            bb.LineWidth = 2;
            pbaspect([1 1 1])
            view([-37 16])
            title(musclenames{ii})
            grid on
            xlabel([jointnames{muscles(2,count2)},' (deg)'])
            ylabel([jointnames{muscles(3,count2)},' (deg)'])
            zlabel('Moment Arm Length (mm)')
            count2 = count2 - 2*mod(ii,2) + 2;
        end
        set(gcf,'Position',[1000 10 800 975])
        saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper\projectedpathsfigure','\threeD','_',datestr(datetime('now'),'yyyymmdd'),'.png'])
        %%Flat graphs
        figure
        count2 = 1;
        for ii = 1:size(muscles,2)
            [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,ii)},axis,muscles(2,ii),0,0);
            jointangles = [-obj.joint_obj{1}.rec_angle_profile(xx)*(180/pi),...
                obj.joint_obj{2}.rec_angle_profile(xx)*(180/pi)-81.578,...
                -obj.joint_obj{3}.rec_angle_profile(xx)*(180/pi)-20.0535];
            momprof = moment_arm_profile(:,3);
            %mnorm = (moment_arm_profile(:,3)-min(moment_arm_profile(:,3)))/(max(moment_arm_profile(:,3))-min(moment_arm_profile(:,3)));
            subplot(3,2,ii)
            if mod(ii,2)==1
                sig = 1;
            else
                sig = -1;
            end
            row = 2;
            bb = plot(sig*jointangles(:,muscles(row,ii)),moment_arm_profile(:,3),'k');
            %bb = fill3(jointangles(:,muscles(3,count2)),jointangles(:,muscles(2,count2)),mnorm,'r');
            hold on
            tdlo = [];
            tdlo(1,1) = sig*jointangles(6,muscles(row,ii));
            tdlo(1,2) = momprof(6);
            tdlo(2,1) = sig*jointangles(51,muscles(row,ii));
            tdlo(2,2) = momprof(51);
            %Touch Down Marker
            td = plot(tdlo(1,1),tdlo(1,2),'ro','MarkerSize',10,'MarkerFaceColor',[1 0 0]);
            %Lift Off Marker
            lo = plot(tdlo(2,1),tdlo(2,2),'go','MarkerSize',10,'MarkerFaceColor',[0 1 0]);
            %cc = [jointangles(:,muscles(2,count2)),jointangles(:,muscles(3,count2)),moment_arm_profile(:,3)];
            cc = [sig*jointangles(:,muscles(row,ii)),moment_arm_profile(:,3)];
            ee = cc(10:20:end,:);
            ff = cc(11:20:end,:);
            C = ff([1;1]*(1:size(ff,1)),:);
            C(1:2:end,:) = ee;
            for j = 1:2:length(C)-1
%                         arrow(C(j,:),C(j+1,:),'Length',5,'Width',0,'TipAngle',30,'Color','b')
                arrow(C(j,:),C(j+1,:),'Length',5,'Width',.001,'Color','b')
                [warnMsg, ~] = lastwarn;
                if ~isempty(warnMsg)
                    arrow FIXLIMITS
                    lastwarn('');
                end
            end
            bb.LineWidth = 2;
            grid on
            xlabel([jointnames{muscles(row,ii)},' (deg)'])
            ylabel('Moment Arm Length (mm)')
            ylim(zlimits(:,ii)')
            pbaspect([1 1 1])
            title(musclenames{ii})
            count2 = count2 - 2*mod(ii,2) + 2;
        end
        set(gcf,'Position',[10 10 800 975])
        saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper\projectedpathsfigure','\flat','_',datestr(datetime('now'),'yyyymmdd'),'.png'])
        keyboard
    end
    if 0
        limits = [9.8 9.8 1.4 1.4 .8 .8
            16 16 7.8 7.8 6 6];
        muscles = [11 11 33 33 27 27;...
            1 2 1 2 2 3];
        musclenames = {'Biceps Femoris Posterior (Hip)';'Biceps Femoris Posterior (Knee)';'Biceps Femoris Anterior (Hip)';...
                            'Biceps Femoris Anterior (Knee)';'Plantaris (Knee)';'Plantaris (Ankle)'};
        for ii = 1:size(muscles,2)
            [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{muscles(1,ii)},axis,muscles(2,ii),0,0);
            jointangles = [-obj.joint_obj{1}.rec_angle_profile(xx)*(180/pi),...
                obj.joint_obj{2}.rec_angle_profile(xx)*(180/pi)-81.578,...
                -obj.joint_obj{3}.rec_angle_profile(xx)*(180/pi)-20.0535];
            %jointangle = moment_arm_profile(:,2);
            momentarm = moment_arm_profile(:,3);
            if muscles(2,ii) ~= 2
                jointname = obj.joint_obj{muscles(2,ii)}.name(4:end-1);
            else
                jointname = obj.joint_obj{muscles(2,ii)}.name(4:end);
            end
            subplot(3,2,ii)
            plot(jointangles(:,muscles(2,ii)),momentarm,'LineWidth',2);
            grid on
            title(musclenames{ii})
            xlabel([jointname,' Joint Angle (deg)'])
            ylabel('Moment Arm Length (mm)')
            pbaspect([1 1 1])
            ylim(limits(:,ii)')
        end
        set(gcf,'Position',[500 20 600 900])
        %saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper','\momentarmloops','_',datestr(datetime('now'),'yyyymmdd'),'.png'])
        keyboard
    end
    for i=1:length(muscles)
%                     figure
%                         subplot(4,1,1)
%                         plot(moment_arm_profile(:,2),'LineWidth',2)
%                         if muscles(1,i) == 33
%                             hold on
%                             plot(obj.theta_motion(xx,2)*(180/pi),'LineWidth',2)
%                         end
%                         grid on
%                         title(obj.musc_obj{muscles(1,i)}.muscle_name(4:end))
%                         ylabel([obj.joint_obj{muscles(2,i)}.name(4:end),' Angle (deg)'])
%                         xlabel('Stride Percent')
%                         subplot(4,1,2)
%                         plot(moment_arm_profile(:,3),'LineWidth',2)
%                         grid on
%                         ylabel('Moment Arm Length (mm)')
%                         xlabel('Stride Percent')
%                         subplot(4,1,3)
%                         plot(moment_arm_profile(:,2),moment_arm_profile(:,3),'LineWidth',2)
%                         grid on
%                         ylabel('Moment Arm Length (mm)')
%                         xlabel([obj.joint_obj{muscles(2,i)}.name(4:end),' Angle (deg)'])
%                         subplot(4,1,4)
%                         anglenorm = (moment_arm_profile(:,2)-min(moment_arm_profile(:,2)))/(max(moment_arm_profile(:,2))-min(moment_arm_profile(:,2)));
%                         momentnorm = (moment_arm_profile(:,3)-min(moment_arm_profile(:,3)))/(max(moment_arm_profile(:,3))-min(moment_arm_profile(:,3)));
%                         plot(anglenorm,'LineWidth',2)
%                         hold on
%                         plot(momentnorm,'LineWidth',2)
%                         grid on
%                         ylabel('Normalized Parameter')
%                         xlabel('Stride Percent')
%                         %legend({'Normalized Angle','Normalized Moment Arm Length'},'Location','southoutside')
%                         set(gcf,'Position',[600 50 500 900])
%                         saveas(gcf,['G:\My Drive\Rat\BiomimeticsPaper','\whole_',obj.musc_obj{muscles(1,i)}.muscle_name(4:end),'_',datestr(datetime('now'),'yyyymmdd'),'.png'])
% %%
%                         for k = 1:3
%                             relevant_muscles = [];
%                             for ii = 1:38
%                                 attachment_bodies = cell2mat(obj.musc_obj{ii}.pos_attachments(:,3));
%                                 if joint == 1 
%                                     if attachment_bodies(1) == 1
%                                         relevant_muscles = [relevant_muscles;ii];
%                                     end
%                                 elseif joint == 2
%                                     if attachment_bodies(end) == 3
%                                         relevant_muscles = [relevant_muscles;ii];
%                                     end
%                                 elseif joint == 3
%                                     if attachment_bodies(end) == 4
%                                         relevant_muscles = [relevant_muscles;ii];
%                                     end
%                                 end
%                             end
%                             figure
%                             for l=1:length(relevant_muscles)
%                                 if l == 1
%                                     [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{relevant_muscles(l)},axis,k,0,0);
%                                     anglenorm = (moment_arm_profile(:,2)-min(moment_arm_profile(:,2)))/(max(moment_arm_profile(:,2))-min(moment_arm_profile(:,2)));
%                                     momentnorm = (moment_arm_profile(:,3)-min(moment_arm_profile(:,3)))/(max(moment_arm_profile(:,3))-min(moment_arm_profile(:,3)));
%                                     deriv = diff(momentnorm);
%                                     if (momentnorm(50)-momentnorm(1)) < 0
%                                         if .5*(momentnorm(40)+momentnorm(60)) > .15
%                                             if momentnorm(1) < .55
%                                                 subplotter = 5;
%                                             else
%                                                 subplotter = 4;
%                                             end
%                                         else
%                                             subplotter = 3;
%                                         end
%                                     else
%                                         subplotter = 2;
%                                     end
%                                     subplot(5,1,1)
%                                     plot(anglenorm,'LineWidth',2)
%                                     subplot(5,1,subplotter)
%                                     plot(momentnorm,'LineWidth',2)
%                                     hold on
%                                 else
%                                     [moment_arm_profile,~] = compute_muscle_moment_arm(obj,obj.musc_obj{relevant_muscles(l)},axis,k,0,0);
%                                     momentnorm = (moment_arm_profile(:,3)-min(moment_arm_profile(:,3)))/(max(moment_arm_profile(:,3))-min(moment_arm_profile(:,3)));
%                                     deriv = diff(momentnorm);
%                                     if (momentnorm(50)-momentnorm(1)) < 0
%                                         if .5*(momentnorm(40)+momentnorm(60)) > .15
%                                             if momentnorm(1) < .55
%                                                 subplotter = 5;
%                                             else
%                                                 subplotter = 4;
%                                             end
%                                         else
%                                             subplotter = 3;
%                                         end
%                                     else
%                                         subplotter = 2;
%                                     end
%                                     subplot(5,1,subplotter)
%                                     plot(momentnorm,'LineWidth',2)
%                                     hold on
%                                 end
%                                 grid on
%                                 ylabel('Normalized Parameter')
%                                 xlabel('Stride Percent')
%                                 %legend({'Normalized Angle','Normalized Moment Arm Length'},'Location','southoutside')
%                             end
%                        end
                %%
    end
end