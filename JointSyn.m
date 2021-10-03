classdef JointSyn < matlab.mixin.SetGet
    properties
        name
        index
        init_pos
        init_pos_w
        init_rot
        Cabs
        CR
        joint_rotmat_profile
        enable_limit = 0;
        euler_angs
        fricCoeff
        limits
        rec_angle_time
        rec_angle_profile
        rec_angleDot_profile
        rec_torque_profile
        sim_position_profile
        type
        uu_joint
        uu_joint2
        uu_joint3
        uuw_joint
        uuw_joint2
        uuw_joint3
    end
    methods
        function obj = JointSyn(Jname)
            obj.name = Jname;
        end
        %% Function: Compute Muscle Force for Motion (outputs Musc_Tension, Musc_Length, Musc_vel, MN_act, etc)
        function [musc_tension, musc_length, musc_velocity, musc_act, MN_act, musc_tension_passive] = compute_musc_force_for_motion(obj,legobj)
            %Using data built into the object, compute the muscle forces
            %necessary to generate known torques at known positions and
            %velocities.
            if isempty(obj.rec_torque_profile) || isempty(obj.rec_angle_profile) || isempty(obj.rec_angleDot_profile)
                disp('call Joint.set_torque_and_kinematic_data()')
            end
            
            %The load on the joint is the opposite of the applied torque
%             torque_load = -obj.rec_angle_profile;
            net_torque = zeros(size(-obj.rec_torque_profile));

            
            %Initialize our outputs
            if isempty(obj.rec_angle_profile)
                T_ext = zeros(obj.num_pts,1);
                T_flx = zeros(obj.num_pts,1);
                obj.L_ext = zeros(obj.num_pts,1);
                obj.L_flx = zeros(obj.num_pts,1);
                num_loops = obj.num_pts;
                obj.rec_angle_profile = obj.theta;
            else
                num_loops = length(obj.rec_angle_profile);
                T_ext = zeros(num_loops,1);
                T_flx = zeros(num_loops,1);
                musc_act = zeros(num_loops,2);
                musc_group_act = zeros(num_loops,2);
                MN_act = -.100 + zeros(num_loops,2);
                MN_group_act = -.1 + zeros(num_loops,2);
                musc_tension_passive = MN_act;
            end
            
            %[muscle_length,muscle_velocity,f_vec] = legobj.musc_obj{1}.compute_muscle_length(legobj);
            
            for i=1:num_loops
                
                %We will now solve for the muscle forces required to
                %counter the maximum torque with the maximum joint
                %stiffness. 
                bnds = [0,max(obj.Amp_ext,obj.Amp_flx)];

                %Explicit solution
                %Flexor passive force
                if i > 1
                    T_passive_flx = (obj.c_flx*T_flx(i-1)/obj.dt_motion/obj.kse_flx + obj.kpe_flx*max(0,flexor_length(i)-obj.l_flx_rest) + obj.c_flx*flexor_velocity(i))/(1+(obj.kpe_flx + obj.c_flx/obj.dt_motion)/obj.kse_flx);
                else
                    T_passive_flx = (obj.kpe_flx*max(0,flexor_length(i)-obj.l_flx_rest) + obj.c_flx*flexor_velocity(i))/(1+(obj.kpe_flx + obj.c_flx)/obj.kse_flx);
                end
                
                ext_passive_muscle_torque = obj.net_joint_torque(T_passive_ext,'ext',obj.rec_angle_profile(i),0,0);
                flx_passive_muscle_torque = obj.net_joint_torque(T_passive_flx,'flx',obj.rec_angle_profile(i),0,0);
                
                net_torque(i) = -obj.rec_torque_profile(i) + ext_passive_muscle_torque + flx_passive_muscle_torque;

                if net_torque(i) > 0
                    %The external load is positive, so the joint must be
                    %supplying a negative torque to cancel it.
                    active_musc_string = 'flx';
                else
                    active_musc_string = 'ext';
                end

                %Now, find the muscle tension that will generate the
                %desired torque by finding the zero of this function.
                f_to_opt = @(x)obj.net_joint_torque(x,active_musc_string,obj.rec_angle_profile(i),net_torque(i),0);
                if isnan(f_to_opt(bnds(1))) || isnan(f_to_opt(bnds(2)))
                    disp('f isnan')
                    keyboard
                end
                T_active = bisect(f_to_opt,bnds(1),bnds(2),1e-6,1e-6,1000);
                
                if net_torque(i) > 0
                    %The flexor is applying the net force
                    T_flx(i) = T_active + T_passive_flx;
                    T_ext(i) = T_passive_ext;
                else
                    %The extensor is applying the net force
                    T_flx(i) = T_passive_flx;
                    T_ext(i) = T_active + T_passive_ext;
                end
                    
                musc_tension_passive(i,:) = [T_passive_ext,T_passive_flx];
                
            end

            musc_tension = [T_ext,T_flx];
            T_ext_dot = [0;diff(T_ext)/obj.dt_motion];
            T_flx_dot = [0;diff(T_flx)/obj.dt_motion];
%             keyboard
            musc_length = [extensor_length,flexor_length];
            musc_velocity = [extensor_velocity,flexor_velocity];
            
            if strcmp(obj.joint_name,'LH_HipZ')
                muscle_num = [10;8];
            elseif strcmp(obj.joint_name,'LH_Knee')
                muscle_num = [3;6];
            elseif strcmp(obj.joint_name,'LH_AnkleZ')
                muscle_num = [9;2];
            end
            
            for i=1:num_loops
                if net_torque(i) > 0  
                    %flexor is providing force
                    
                    %Solve for the activation of the muscle by rearranging
                    %dT/dt
                    Act_flx = max(0,((obj.c_flx/obj.kse_flx)*T_flx_dot(i) + (obj.kse_flx+obj.kpe_flx)*T_flx(i))/obj.kse_flx - obj.kpe_flx*max(0,(flexor_length(i)-obj.l_flx_rest)) - obj.c_flx*flexor_velocity(i));  
%                   Act_flx = max(0,(obj.c_flx*T_flx_dot(i) + (obj.kse_flx+obj.kpe_flx)*T_flx(i))/obj.kse_flx - obj.kpe_flx*max(0,(flexor_length(i)-obj.l_flx_rest)) - obj.c_flx*flexor_velocity(i));
%                   Act_flx = max(0,(obj.c_flx*T_flx_dot(i) + (obj.kse_flx+obj.kpe_flx)*T_flx(i))/obj.kse_flx - obj.kpe_flx*(flexor_length(i)-obj.l_flx_rest) - obj.c_flx*flexor_velocity(i)); 
                    musc_act(i,2) = Act_flx;
                    musc_group_act(i,2) = Act_flx/muscle_num(2);
                    flx_sigmoid = @(V) (obj.Amp_flx./(1+exp(obj.steep_flx*(obj.xoff_flx-V))) + obj.yoff_flx)*max(0,(1-(flexor_length(i)-obj.l_flx_rest)^2/(obj.l_width_flx)^2));
                    f_act_flx = @(V) Act_flx - (0*(V <= -.100) + obj.Amp_flx*(V >= 0) + (V>-.100 && V < 0)*flx_sigmoid(V));
                    MN_act(i,2) = bisect(f_act_flx,-.100,0,1e-6,1e-6,1000);
                    MN_group_act(i,2) = ((MN_act(i,2)+.1)/muscle_num(2))-.1;
                else
                    Act_ext = max(0,((obj.c_ext/obj.kse_ext)*T_ext_dot(i) + (obj.kse_ext+obj.kpe_ext)*T_ext(i))/obj.kse_ext - obj.kpe_ext*max(0,(extensor_length(i)-obj.l_ext_rest)) - obj.c_ext*extensor_velocity(i));
%                   Act_ext = max(0,(obj.c_ext*T_ext_dot(i) + (obj.kse_ext+obj.kpe_ext)*T_ext(i))/obj.kse_ext - obj.kpe_ext*max(0,(extensor_length(i)-obj.l_ext_rest)) - obj.c_ext*extensor_velocity(i));
%                   Act_ext = max(0,(obj.c_ext*T_ext_dot(i) + (obj.kse_ext+obj.kpe_ext)*T_ext(i))/obj.kse_ext - obj.kpe_ext*(extensor_length(i)-obj.l_ext_rest) - obj.c_ext*extensor_velocity(i));
                    musc_act(i,1) = Act_ext;
                    musc_group_act(i,1) = Act_ext/muscle_num(1);
                    %Calculate the tension generated by a stimulus. This is
                    %the Animatlab "Stimulus-Tension" curve
                    ext_sigmoid = @(V) (obj.Amp_ext./(1+exp(obj.steep_ext*(obj.xoff_ext-V))) + obj.yoff_ext)*max(0,(1-(extensor_length(i)-obj.l_ext_rest)^2/(obj.l_width_ext)^2));
                    f_act_ext = @(V) Act_ext - (0*(V <= -.100) + obj.Amp_ext*(V >= 0) + (V>-.100 && V < 0)*ext_sigmoid(V));
                    MN_act(i,1) = bisect(f_act_ext,-.100,0,1e-6,1e-6,1000);
                    MN_group_act(i,1) = ((MN_act(i,1)+.1)/muscle_num(1))-.1;
                end
                
            end
            %%
            if 0
                xx = linspace(0,100,num_loops/3)';
                width = 8;
                fontsize = 50;
                fontsizeax = 45;
                figW = 2000;
                %figH = figW/1.75;
                figH = 1050;
                lshift = 10;
                fig1 = figure('pos',[0 0 figW figH]);
                title(strcat(strrep(obj.joint_name,'LH_',' '),' MN Group Activation'),'fontsize',fontsize)
                %view([-38.4 42]);
                ylabel('Percent Stride (%)','fontsize',fontsize)
                view([90 0]);
                %grid off
                hold on
                for ii = 1:muscle_num(1)
                    zz = ii*ones(size(xx));
    %                qq = area(xx,ii*musc_group_act(:,1),'FaceColor',[0 0 0]);
    %                plot3(zz,xx,qq)
                    %plot3(zz,xx,ii*musc_group_act(:,1),'b');
                    %%Exponential
                    %plot3(zz,xx,(ii*(MN_group_act(:,1)+.1)-.1),'Color',[1-(ii*1/muscle_num(1))^2 1-(ii*1/muscle_num(1))^2 1])
                    %%Linear
                    kk = plot3(zz,xx,1000*(ii*(MN_group_act(1:100,1)+.1)-.1),'Color',[1-ii*1/muscle_num(1) 1-ii*1/muscle_num(1) 1]);
                end

                for ii = 1:muscle_num(2)
                    zz = ii*ones(size(xx));
                    %plot3(zz,xx,ii*musc_group_act(:,2),'r')
                    %%Exponential
                    %plot3(zz,xx,(ii*(MN_group_act(:,2)+.1)-.1),'Color',[1 1-(ii*1/muscle_num(2))^2 1-(ii*1/muscle_num(2))^2])
                    kk = plot3(zz,xx,1000*(ii*(MN_group_act(1:100,2)+.1)-.1),'Color',[1 1-ii*1/muscle_num(2) 1-ii*1/muscle_num(2)]);
                end
                    set(findall(gca, 'Type', 'Line'),'LineWidth',width);
                    set(gca,'FontSize',fontsizeax)
                    zlh = zlabel('MN Activation (mV)','fontsize',fontsize);
%                    set(get(gca,'zlabel'),'VerticalAlignment','middle')
%                     zlabpos = get(zlh,'Position');
%                     zlabpos(2) = zlabpos(2)-lshift;
%                     set(zlh,'Position',zlabpos);

                figure('pos',[0 0 figW figH])
                title(strcat(strrep(obj.joint_name,'LH_',' '),' Muscle Group Tension'),'fontsize',fontsize)
                %view([-38.4 42]);
                ylabel('Percent Stride (%)','fontsize',fontsize)
                view([90 0]);
                %grid off
                hold on
                for ii = 1:muscle_num(1)
                    zz = ii*ones(size(xx));
    %                 qq = area(xx,ii*musc_group_act(:,1),'FaceColor',[0 0 0]);
    %                 plot3(zz,xx,qq)
                    %plot3(zz,xx,ii*musc_group_act(:,1),'b');
                    %%Exponential
                    %plot3(zz,xx,(ii*musc_group_act(:,1)),'Color',[1-(ii*1/muscle_num(1))^2 1-(ii*1/muscle_num(1))^2 1])
                    %%Linear
                    pp = plot3(zz,xx,(ii*musc_group_act(1:100,1)),'Color',[1-ii*1/muscle_num(1) 1-ii*1/muscle_num(1) 1]);

                end

                for ii = 1:muscle_num(2)
                    zz = ii*ones(size(xx));
                    %plot3(zz,xx,ii*musc_group_act(:,2),'r')
                    %%Exponential
                    %plot3(zz,xx,(ii*musc_group_act(:,2)),'Color',[1 1-(ii*1/muscle_num(2))^2 1-(ii*1/muscle_num(2))^2])
                    %%Linear
                    pp = plot3(zz,xx,(ii*musc_group_act(1:100,2)),'Color',[1 1-ii*1/muscle_num(2) 1-ii*1/muscle_num(2)]);
                end
                    zlh = zlabel('Muscle Activation (N)','fontsize',fontsize);
                    set(findall(gca, 'Type', 'Line'),'LineWidth',width);
                    set(gca,'FontSize',fontsizeax)
%                    set(get(gca,'zlabel'),'VerticalAlignment','middle')
%                     zlabpos = get(zlh,'Position');
%                     zlabpos(2) = zlabpos(2)-lshift;
%                     set(zlh,'Position',zlabpos);
                
                uu = figure('pos',[0 0 figW figH]);
%                 axes1 = axes('Parent',uu,...
%                     'Position',[0.251423921887714 0.231477657085432 0.653576078112285 0.638281796312414]);
                hold on
                for ii = 1:size(MN_act,2)
                    ooo = plot(xx,1000*MN_act(1:100,ii));
                    if ii == 1
                        set(ooo,'Color',[0 0 1]);
                    else
                        set(ooo,'Color',[1 0 0]);
                    end
                end
                title('Hip Generalized MN Activation','fontsize',fontsize);
                xlabel('Percent Stride (%)','fontsize',fontsize);
                ylh = ylabel('MN Activation (mV)','fontsize',fontsize);
                set(gca,'FontSize',fontsizeax);
%                set(get(gca,'ylabel'),'VerticalAlignment','middle');
%                 ylabpos = get(ylh,'Position');
%                 ylabpos(1) = ylabpos(1)-lshift;
%                 set(ylh,'Position',ylabpos);
                legend({'Extension','Flexion'},'FontSize',fontsize);
                set(findall(gca, 'Type', 'Line'),'LineWidth',width);
                
                ll = figure('pos',[0 0 figW figH]);
%                 axes1 = axes('Parent',ll,...
%                     'Position',[0.209260908281389 0.231477657085432 0.695739091718609 0.638281796312414]);
                hold on
                for ii = 1:size(musc_act,2)
                    bbb = plot(xx,musc_act(1:100,ii));
                    if ii == 1
                        set(bbb,'Color',[0 0 1]);
                    else
                        set(bbb,'Color',[1 0 0]);
                    end
                end
                title('Hip Generalized Muscle Activation','fontsize',fontsize);
                xlabel('Percent Stride (%)','fontsize',fontsize);
                ylh = ylabel('Muscle Activation (N)','fontsize',fontsize);
                set(gca,'FontSize',fontsizeax);
%                set(get(gca,'ylabel'),'VerticalAlignment','middle')
%                 ylabpos = get(ylh,'Position');
%                 ylabpos(1) = ylabpos(1)-lshift;
%                 set(ylh,'Position',ylabpos);
                legend({'Extension','Flexion'},'FontSize',fontsize);
                set(findall(gca, 'Type', 'Line'),'LineWidth',width);
                
                %fname = 'G:\My Drive\Rat\Optimizer\OutputFiguresGroup';
                fname = 'G:\My Drive\Rat\Optimizer\OutputFiguresGroup';
                saveas(kk, fullfile(fname,strcat(obj.joint_name,' MN Group Activation')), 'png');
                saveas(pp, fullfile(fname,strcat(obj.joint_name,' Muscle Group Tension')), 'png');
                saveas(uu, fullfile(fname,strcat(obj.joint_name,' Generalized MN Activation')), 'png');
                saveas(ll, fullfile(fname,strcat(obj.joint_name,' Generalized Muscle Activation')), 'png');
                close all
            end
                %keyboard
            %%
            %Now that we have muscle lengths, velocities, and forces for
            %every time step in our data, we can calculate what our
            %activation and passive parameters should look like
            for i=1:obj.num_pts
                %We will now solve for the muscle forces required to
                %counter the maximum torque with the maximum joint
                %stiffness. This is a smooth function, so we will use our
                %quasi-newton optimizer. We will constrain it such that
                %muscle forces are greater than 0. A large upper bound will
                %be used since none is not an option. 
                bnds = [0,1e6];

                %First, assume the torque is positive. This will reveal the
                %maximum flexor force.
                f_to_opt = @(x)obj.net_joint_torque(x,'ext',obj.theta(i),obj.rec_torque_profile(i),1);
                [T,ext_res] = con_opti(f_to_opt,[1],{[],[]},bnds,[],[],[1000,1e-12,1e-12],'min','bfgs',[],'linesearch',{'backtrack',[1e-4,.9,.5,500],[],[]},[],0);
                T_ext_ext(i) = T(1);
                T_flx_ext(i) = 0;
                
                obj.L_ext(i) = obj.l_ext_current;
                
                f_to_opt = @(x)obj.net_joint_torque(x,'flx',obj.theta(i),-obj.rec_torque_profile(i),1);
                [T,flx_res] = con_opti(f_to_opt,[1],{[],[]},bnds,[],[],[1000,1e-12,1e-12],'min','bfgs',[],'linesearch',{'backtrack',[1e-4,.9,.5,500],[],[]},[],0);
                T_ext_flx(i) = 0;
                T_flx_flx(i) = T(1);

                obj.L_flx(i) = obj.l_flx_current;
                
                force_residuals(i,:) = [ext_res,flx_res];
            end
            
            %Stores the minimum lengths calculated into the properties
            obj.min_l_ext = min(obj.L_ext);
            obj.max_l_ext = max(obj.L_ext);
            obj.min_l_flx = min(obj.L_flx);
            obj.max_l_flx = max(obj.L_flx);
            
            obj.MN_ext = MN_act(:,1);
            obj.MN_flx = MN_act(:,2);
        end
    end
end