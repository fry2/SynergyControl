function [Am_musc,V_musc,Am_musc_raw] = Am_generator(obj,forces)
    % For input forces, generate Am and Vm values for all muscles
    % Input: forces: array of muscle forces over time
    % Output: Am_musc: array of Am values over time
    % Output: V_musc: array of membrane voltages over time

    if size(forces,2)<size(forces,1)
        forces = forces';
    end
%     [beg,ennd,~] = find_step_indices(obj);
    beg = obj.sampling_vector(1);
    ennd = obj.sampling_vector(end);
    odt = obj.dt_motion;
    dt = ((ennd-beg)*odt)/length(forces);
    [forces_dot,forces_doty] = gradient(forces,dt);
    numMuscles = length(obj.musc_obj);

    fl = @(Lm,Lr,Lw) max(1-((Lm-Lr).^2./Lw^2),.7);
%     Am = @(Al,b,ks,T_dot,kp,delL,L_dot,T,Fmax) (1./Al).*((b./ks).*T_dot+(1+kp./ks).*T-min(max(kp.*delL+b.*L_dot,0),Fmax));
    Am = @(Al,b,ks,T_dot,kp,delL,L_dot,T,Fmax) (1./Al).*((b./ks).*T_dot+(1+kp./ks).*T-kp.*delL-b.*L_dot);
    %Am = @(Al,b,ks,T_dot,kp,T,pt) (1./Al).*((b./ks).*T_dot+(1+kp./ks).*T-pt');
    %force_eq = @(Al,b,ks,Am,kp,delL,L_dot,T) (ks./b).*(kp.*delL+b.*L_dot-(1+kp./ks).*T+Al.*Am);
    %Am = @(Al,b,ks,T_dot,kp,delL,L_dot,T) (1./Al).*((b./ks).*T_dot+(1+kp./ks).*T);
    V = @(xoff,ST_max,steepness,yoff,Am) xoff-(1./steepness).*log(((ST_max)./(Am-yoff))-1);
    st_curve = @(ST_max,steepness,xoff,V,yoff) (ST_max./(1+exp(steepness*(xoff-V))))+yoff;
    ssTeqn = @(Am,ks,kp,mL,Lr,Al) (ks/(ks+kp)).*max(kp.*(mL-Lr)+Al.*Am,0);
    Am_musc = zeros(size(forces));
    Am_musc_raw = Am_musc;
    Al_musc_all = Am_musc;
    V_musc = Am_musc;

    for ii = 1:numMuscles
       [b,ks,kp,Lw,Lr,xoff,Fmax,steepness,mL,mV,ST_max,yoff] = getMuscParams(obj,ii,beg,ennd);
       delL_musc = max(mL-Lr,0);
       Al_musc = fl(mL,Lr,Lw);
       %Tdot = forces_dot(ii,:);
       Tdot = smoothdata(forces_dot(ii,:),'movmedian',50);
       T = forces(ii,:);
       Am_musc_raw(ii,:) = Am(Al_musc',b,ks,Tdot,kp,delL_musc',mV',T,Fmax);
       %Am_musc_raw(ii,:) = Am(Al_musc',b,ks,Tdot,kp,T,obj.passive_tension(:,ii));
       %passiveForces = kp*delL_musc+b*mV; passiveForces(passiveForces < 0) = 0; passiveForces(passiveForces > Fmax) = Fmax; 
       Am_musc(ii,:) = min(max(Am_musc_raw(ii,:),0),ST_max+yoff);
       Am_musc(ii,:) = smoothdata(Am_musc(ii,:),'movmedian',50);
       %Am_musc(ii,:) = smoothdata(Am_musc(ii,:),'gaussian',150);
       V_musc(ii,:) = real(V(xoff,ST_max,steepness,yoff,Am_musc(ii,:)));
            V_musc(ii,isinf(V_musc(ii,:))) = -.04;
       V_musc(ii,V_musc(ii,:)<-.06) = -.06;
       V_musc(ii,V_musc(ii,:)>-.04) = -.04;
       Al_musc_all(ii,:) = Al_musc;
       %% For loop plotter 1: Plot 5 subplot fig of tension equation
       if 0
       figure
        subplot(5,1,1)
            plot((1+kp./ks).*T(1:end-10))
            xlabel('(1+kp./ks).*T')
        subplot(5,1,2)
            plot((b./ks).*Tdot(1:end-10))
            xlabel('(b./ks).*Tdot')
        subplot(5,1,3)
            plot(kp.*delL_musc(1:end-10))
            xlabel('kp.*delL_musc','Interpreter','None')
        subplot(5,1,4)
            plot(b.*mV(1:end-10))
            xlabel('b.*mV')
        subplot(5,1,5)
            plot(Am_musc(ii,1:end-10),'LineWidth',2)
            hold on
            plot(zeros(length(Am_musc(ii,1:end-10)),1))
            xlabel('Am_musc','Interpreter','None')
            %ylim([0 max(Am_musc(ii,:))*1.1])
       end
       %% For loop plotter 2: View one muscle's active and passive waveforms. When active is > passive, Am is possible.
       if 0
            active = (b/ks).*Tdot+(1+kp/ks).*T;
            passive = kp.*delL_musc+b.*mV;
            %yLims = [min([active',passive],[],'all') max([active',passive],[],'all')];
            figure
            plot(active(1:end-10),'r','LineWidth',2)
            title('Am possible where active >= passive')
            hold on
            plot(passive(1:end-10),'b','LineWidth',2)
            legend({'Active','Passive'})
            %ylim(yLims)
       end
       %%
    end
    
    % Low values of necessary Am output cause rapid shifts in applied voltage. Rather than simply cap the voltages at a bottom value and
    ... creating many flat plateaus, we smooth the data. This is preferable over a moving lowpass filter because it doesn't distort the ends of the data
    %V_musc = smoothdata(V_musc','gaussian',10)';

    if 0
        figure
        plot(outT')
        title('Setting Am < 0 to zero')
        figure
        plot(forces')
        title('Original forces trying to recreate')
    end
end
    %% inputAmoutputT
    function outT = inputAmoutputT(obj,mnum,mL,mV,Am)
        [b,ks,kp,Lw,Lr,xoff,Fmax,steepness,mL,mV,ST_max,yoff] = getMuscParams(obj,mnum,1,length(obj.theta_motion));
        dt = obj.dt_motion;
        Tdot(1) = (ks/b)*(kp*max(mL(1)-Lr,0));
        iT(1) =  Tdot(1)*dt;
        T(1) = 0;
        Al = @(Lm,Lr,Lw) max(1-((Lm-Lr).^2./Lw^2),.7);
        
        for i = 2:length(mV)
            Tdot(i) =  (ks/b)*(kp*max(mL(i)-Lr,0) + b*mV(i) - (1+kp/ks)*iT(i-1) + Al(mL(i),Lr,Lw)*Am(i));  

            %The new tension
            iT(i) =  iT(i-1) +  Tdot(i)*dt;

            if iT(i) >= 0
                T(i) =  iT(i);
            else
                T(i) = 0;
            end

            %Make certain that the tension never exceed the absolute maximum set by the user.

            if( T(i) >  Fmax)
                T(i) =  Fmax;
            end
        end

        outT = T;
    end
    %% getMuscParams
    function [b,ks,kp,Lw,Lr,xoff,Fmax,steepness,mL,mV,ST_max,yoff] = getMuscParams(obj,mnum,beg,ennd)
       m = length(obj.theta_motion(beg:ennd));
       n = length(obj.sampling_vector);
       musc = obj.musc_obj{mnum};
       b = musc.damping;
       ks = musc.Kse;
       kp = musc.Kpe;
       Lw = musc.l_width;
       Lr = musc.RestingLength;
       xoff = musc.x_off;
       Fmax = musc.max_force;
       steepness = musc.steepness;
       mL = musc.muscle_length_profile;
       mV = musc.muscle_velocity_profile;
       ST_max = musc.ST_max;
       yoff = musc.y_off;
    end
    %% findParamsfromLloyd
    function x = findParamsFromLloyd(Tdot,Al_musc,T,mV,mL,Lr)
        % Try to find ks, b, kp, Am values that satisfy an input force waveform
        % Doesn't really work atm
        lloydopt = @(x) lloyd2animatlab(x,Tdot',Al_musc,T',mV,mL,Lr);

        % Find possible muscle parameter solutions
        rng default % for reproducibility
        N = 100; % try 100 random start points
        pts = 1000*rand(N,4);
        soln = zeros(N,4); % allocate solution
        fval = zeros(N,1); % allocate solution
        opts = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','MaxFunctionEvaluations',10e6,'PlotFcn',@optimplotfirstorderopt,...
            'FunctionTolerance',1);

        lb = [0,0,0,0];
        ub = [10e4,10e4,10e4,10e4];
        A = [];
        B = [];
        Aeq = [];
        beq = [];
        x0 = [1000,10,100,10];
        mL_rel = max(mL-Lr,0);
        lloydopt = @(x) lloyd2animatlab(x,Tdot(1),Al_musc(1),T(1),mV(1),mL(1),Lr);
        fun = @(x) Tdot(1)-(x(1)./x(2)).*(x(3).*mL_rel(1)+x(2).*mV(1)-(1+x(3)./x(1)).*T(1)+x(4).*Al_musc(1));
        [x,fval,exitflag,output] = fmincon(lloydopt,x0,A,B,Aeq,beq,lb,ub);

        for k = 1:N
            [soln(k,:),fval(k,1)] = fsolve(lloydopt,pts(k,:),opts); % find solutions
        end
    end