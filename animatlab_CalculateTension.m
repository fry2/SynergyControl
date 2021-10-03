function [Tension,Tdot,Act,mL,h] = animatlab_CalculateTension(obj,muscle,Kpe,to_plot) 
    mL = muscle.muscle_length_profile;
    %simLens = importdata([pwd,'\Animatlab\SynergyWalking\DataTool_1.txt']);
    %mL = simLens.data(:,contains(simLens.colheaders,'mL'));
    Lr = muscle.RestingLength;
    amp = (.5*Lr);
    offset = Lr;
    %mL = amp*sin(3.2.*linspace(0,2*pi,length(mL)))+offset;
    %mL = linspace(1.5*Lr,.5*Lr,length(mL));
    %mL = .7*Lr.*ones(1,length(mL));
    dt = obj.dt_motion;
    
    %temp = importdata([pwd,'\Animatlab\SynergyWalking\DataTool_2.txt']);
    %ndat = temp.data;
    
    m_iVelAvgIndex = 1;
    m_iMuscleVelAvgCount = 4;
    
    % Constants
    PELenPerc = .9;
    MinPeLenPerc = .005;
    SeRestLength = Lr - (Lr*PELenPerc);
    PeMin = Lr*MinPeLenPerc;
    Vm = -.050;
    %Iapp = linspace(20,0,length(mL));
    Iapp = zeros(1,length(mL));
    Iapp(3705:4631) = 10;
    time = 0:dt:(length(mL)-1)*dt;
    tau = .7e-3; %s
    taums = .7;
    dtms = dt*1e3;
    MaxTension = muscle.max_force;
    %Kpe = muscle.Kpe;
    Kse = muscle.Kse;
    
    B = muscle.damping;
    KseByB = Kse/B;
    KpeByKse = (1 + (Kpe/Kse));
    
    % Initializations
    Vse(1) = 0;
    Vpe(1) = 0;
    Vmuscle(1:2) = 0;
    Displacement(1) =  mL(1)-Lr;
    Tdot(1) = KseByB*( Kpe* Displacement(1));
    InternalTension(1) =  Tdot(1)*dt;
    Tension(1) = InternalTension(1);
    SeLength(1) = SeRestLength + Tension(1)/Kse;
    PeLength(1) = mL(1)-SeLength(1);
    U(1) = 1e-3*exp((-dtms*taums)/(1-dtms))*Iapp(1);
    Usim = @(t) max(Iapp)*(1-exp(-t/tau));
    U(1) = Usim(0);
    Vm(1) = U(1)-.06;
    TL(1) = Ftl(mL(1),Lr);
    
    for i = 2:length(mL)
        %Store the previous muscle length
         PrevLength = mL(i-1);

        %Calculate the current muscle length.
         Length(i) = mL(i);

        %Calculate the displacement of this muscle d = (x-x*)
         Displacement(i) =  Length(i)-Lr;
         DisplacementRatio(i) =  Length(i)/Lr;

        %Calculate the instantaneous velocity of change of the muscle length.
        if i > 2
%            Vmuscle(i) = (mL(i-1)- mL(i-2))/dt;
             Vmuscle(i) = (Length(i)- PrevLength)/dt;
        end

        %Calculate the active force that is generated by neural stimulation
         TL(i) = Ftl(Length(i),Lr);
         U(i) = 1e-3*((1-dtms/taums)*U(i-1)+(dtms/taums)*Iapp(i));
         Vm2(i) = (U(i)-.06);
         Vm(i) = ((Iapp(i)/1000)-.06);
         Act(i) = Fact(Vm(i),muscle);
         A(i) =  TL(i)*Act(i);

         TLPerc =  TL(i)*100;
         
        %Calculation of the derivitave of the tension
        
         Tdot(i) =  KseByB*( max(Kpe* Displacement(i),0) +  B*Vmuscle(i) -  KpeByKse*InternalTension(i-1) +  A(i));  

        %The new tension
            InternalTension(i) =  InternalTension(i-1) +  Tdot(i)*dt;

        %tension can never be negative, but we want to maintain the "internal" calculations so that the
        %time constants are correct. If you shorten the muscle rapidly it will take it some time to 
        %rebuild its tension. This time will be determined by the internal tension, new length,  and the params like viscosity.
        %But the force seen on the muscle itself will still be 0 N because it is slack.
        %If we did not do this then if your muscle is at rest and is pulled, stays steady, and then relaxed, then
        %it would end up generating a negative tension. (See fig 4 of shadmehr web doc on muscle model).
        %There is a problem in this figure in that you can see that the muscle produced negative tension and this is not 
        %possible in a real muscle.
        %Tension(i) =  InternalTension(i);
        %Tension(i) =  InternalTension(i);
        if( InternalTension(i) >= 0)
            Tension(i) =  InternalTension(i);
        else
            Tension(i) = 0;
        end

        %Make certain that the tension never exceed the absolute maximum set by the user.
         
        if( Tension(i) >  MaxTension)
            Tension(i) =  MaxTension;
        end

        SeLPrev =  SeLength(i-1);
        PeLPrev = Lr - SeLPrev;

        SeLength(i) =  SeRestLength + (Tension(i)/Kse);
        SeDisplacement(i) =  SeLength(i) -  SeRestLength;
        if( SeDisplacement < 0) 
            SeDisplacement = 0;
        end

        PeLength(i) =  Length(i) -  SeLength(i);
        if( PeLength(i) < PeMin)
            SeLength(i) =  Length(i)  - PeMin;
            PeLength(i) = PeMin;
        end

         Vse(i) = ( SeLength(i)- SeLPrev)/dt;
         Vpe(i) = ( PeLength(i)- PeLPrev)/dt;
    end

    mLr = Lr.*ones(1,length(mL));
    mMax = 1.5*Lr.*ones(1,length(mL));
    mMin = .5*Lr.*ones(1,length(mL));
    
    kp = Kpe;
    a = Act./(Lr*kp);
    badlineL = Lr.*(1+(1-sqrt(1+16.*a.^2))./(8.*a));
    if ~all(size(mL)==size(TL))
        mL = mL';
    end
    badlineAm = (kp.*(Lr-mL))./(TL);
    badlineAm(badlineAm>muscle.ST_max) = muscle.ST_max;
    badlineAm(badlineAm<0) = 0;
    time = 0:dt:(length(mL)-1)*dt;
    if max(time) ~= 10
        simTime = max(time);
    else
        simTime = 10;
    end
    
    fAlpha = .4;
    Act(1) = Act(2);
    titleSize = 16;

%% Big Fig 1    
if to_plot
    figure('Position',[27,34,1195,750]);
    plot1 = subplot(3,1,1);
        h(1) = plot(Iapp,(mL./Lr),'LineWidth',2,'Color','k');
        plot1.FontSize = 12;
        plot1.YTickLabel = {'.5','L_{rest}','1.5'};
        plot1.FontWeight = 'Bold';
        %ylabel('^{L}/_{L_{rest}}','FontSize',14)
        plot1.YLabel.Rotation = 0;
        p1ypos = plot1.YLabel.Position;
        p1ypos(1) = p1ypos(1)-1.2;
        plot1.YLabel.Position = p1ypos;
        xlabel("Applied Current (nA)"+newline+" ",'FontSize',12)
        %xlim([0 max(Iapp)])
        ylim([.5 1.6])
        hold on
        h(2) = area(Iapp,badlineL/Lr,'FaceColor','r','FaceAlpha',fAlpha,'EdgeAlpha',0);
        h(3) = plot(Iapp,badlineL/Lr,'Color','r','LineWidth',2);
        legend(h([1 2]),{'Muscle Length';...
                [newline 'Discontinuous' newline 'Tension' newline 'Region']},...
                'Location','eastoutside',...
                'FontSize',12)
        %title(['Muscle Length ',num2str(Kpe/obj.musc_obj{2}.Kpe)],'FontSize',titleSize)
        drawnow
        pos1 = get(plot1,'Position');
    plot3 = subplot(3,1,2);
        h(1) = plot(Iapp,Act/max(Act),'LineWidth',2,'Color','k');
        plot3.FontSize = 12;
        plot3.YTickLabel = {'0','.5','ST_{max}'};
        plot3.FontWeight = 'Bold';
        %ylabel('^{A_{m}}/_{ST_{max}}','FontSize',14)
        plot3.YLabel.Rotation = 0;
        p3ypos = plot3.YLabel.Position;
        p3ypos(1) = p3ypos(1)-1.2;
        plot3.YLabel.Position = p3ypos;
        xlabel("Applied Current (nA)"+newline+" ",'FontSize',12)
        %xlim([0 max(Iapp)])
        hold on
        h(2) = area(Iapp,badlineAm/max(Act),'FaceColor','r','FaceAlpha',fAlpha,'EdgeAlpha',0);
        h(3) = plot(Iapp,badlineAm/max(Act),'Color','r','LineWidth',2);
        legend(h([1 2]),{'Activation';...
                [newline 'Discontinuous' newline 'Tension' newline 'Region']},...
                'Location','eastoutside',...
                'FontSize',12);
        pos3 = plot3.Position;
        %set(plot3,'Position',[pos3(1:2),pos1(3:4)])
        title('Activation','FontSize',titleSize)
        drawnow
    plot4 = subplot(3,1,3);
        plot(Iapp,Tension,'LineWidth',2,'Color','k')
        plot4.FontSize = 12;
        plot4.FontWeight = 'Bold';
        %ylabel('^{Tension}/_{F_{max}}','FontSize',14)
        plot4.YLabel.Rotation = 0;
        p4ypos = plot4.YLabel.Position;
        p4ypos(1) = p4ypos(1)-1.5;
        p4ypos(2) = p4ypos(2)-2;
        plot4.YLabel.Position = p4ypos;
        xlabel('Applied Current (nA)','FontSize',12)
        xlim([0 max(Iapp)])
        hold on 
        tenshade = max(Tension).*double(Tension==0);
        area(Iapp,tenshade,'FaceColor','r','FaceAlpha',fAlpha,'EdgeAlpha',0)
        legend({'Tension';...
                [newline 'Discontinuous' newline 'Tension' newline 'Region']},...
                    'Location','eastoutside',...
                    'FontSize',12)
        ylim([0 max(Tension)])
        pos4 = plot4.Position;
        plot4.YTick = [0,.5*max(Tension),max(Tension)];
        plot4.YTickLabel = {'0','.5','F_{max}'};
        set(plot4,'Position',[pos4(1:2),pos1(3:4)])
        title('Tension','FontSize',titleSize)
        drawnow
end
%%
        
    T_out = [];
    T_find = [];
    T = [];
%    [T_find,T_out,T] = compute_muscle_passive_tension(obj,muscle,0);

	 %IbRate =  IbDischargeConstant* Tension;

%     if ~m_bEnabled
%         Tension = 0;
%         InternalTension = 0;
%     end

    function fltTl = Ftl(fltL, RestingLength)
        fltTl = CalculateLTGain(fltL, RestingLength);
        if fltTl<0
            fltTl = 0;
        end
    end

    function fltTl = CalculateLTGain(fltInput, RestingLength)
        fltLceNorm = fltInput -  RestingLength;
         TLc = .5* RestingLength;
        fltTl = (-((fltLceNorm^2)/ TLc^2)  + 1);
        if fltTl<0 
            fltTl = 0;
        end
    end

    function fltAct = Fact(fltStim,muscle)
        fltAct = CalculateSTGain(fltStim,muscle);
        if fltAct < 0
            fltAct = 0;
        end
    end

    function outVal = CalculateSTGain(fltInput,muscle)
         Bb = muscle.ST_max;
         Aa = muscle.x_off;
         C = muscle.steepness;
         D = muscle.y_off;
        if(InLimits(fltInput,'ST',muscle))
            outVal = (( Bb/(1+exp( C*( Aa-fltInput)))) +  D);
        else
            outVal = CalculateLimitOutput(fltInput,'ST',muscle);
        end
    end

    function outVal = InLimits(fltInput,curveType,muscle)
        m_bUseLimits = 1;
        switch curveType
            case 'LT'
                 LowerLimit = min(muscle.muscle_length_profile);
                 UpperLimit = max(muscle.muscle_length_profile);
            case 'ST'
                 LowerLimit = -.06;
                 UpperLimit = -.04;
        end
        if( m_bUseLimits && ( (fltInput <  LowerLimit) || (fltInput >  UpperLimit) ) ) 
            outVal = false;
        else
            outVal = true;
        end
    end

    function outVal = CalculateLimitOutput(fltInput,curveType,muscle)
        switch curveType
            case 'LT'
                 LowerLimit = min(muscle.muscle_length_profile);
                 UpperLimit = max(muscle.muscle_length_profile);
            case 'ST'
                 LowerLimit = -60;
                 UpperLimit = -40;
        end
        if(fltInput <  LowerLimit)
            outVal =  LowerOutput;
        elseif (fltInput >  UpperLimit)
            outVal =  UpperOutput;
        else
            outVal = 0;
        end
    end
end