function sim_analysis(muscnum,obj,sdata,V_musc,forces)
    variables = whos;
    varNames = {variables(:).name};
    scrSz = get(groot, 'ScreenSize');
    scW = scrSz(3);
    close all

    if ~any(contains(varNames,'sdata'))
        sdata = processSimData("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim",1);
    end

    terms2find = {'SynergyStims';'KeyMNs';'KeyMuscles';'KeyMuscTen';'KeyMuscAct'};
    for ii=1:length(terms2find)
        rowNames = {sdata(:).name};
        try
            rowInd(ii) = find(contains(rowNames,terms2find{ii}));
        catch
            rowInd(ii) = 0;
        end
    end
    clear terms2find rowNames

    syns = sdata(1).data(85:end-10,:);
    keymns = sdata(1).data(85:end-10,:);
    keymuscs = sdata(2).data(85:end-10,:);
    keytens = sdata(3).data(85:end-10,:);
    keyact = sdata(4).data(85:end-10,:);
    n = length(keytens);
    m = 500;

    keymns = interp1(1:length(keymns),keymns,linspace(1,length(keymns),m));
    keymuscs = interp1(1:length(keymuscs),keymuscs,linspace(1,length(keymuscs),m));
    keytens = interp1(1:length(keytens),keytens,linspace(1,length(keytens),m));
    keyact = interp1(1:length(keyact),keyact,linspace(1,length(keyact),m));

    %musc = 2;
    mnInd = find(ismember(sdata(1).data_cols,['neur',num2str(muscnum)]));
    mcInd = find(contains(sdata(2).data_cols,obj.musc_obj{muscnum}.muscle_name(4:end)));
    tnInd = find(contains(sdata(3).data_cols,obj.musc_obj{muscnum}.muscle_name(4:end)));
    acInd = find(contains(sdata(4).data_cols,obj.musc_obj{muscnum}.muscle_name(4:end)));

    xax = linspace(0,100,500);

    [ha, pos] = tight_subplot(3, 1, [.03 .03],[.05 .05],[.08 .05]);

    % subplot(4,2,1)
    % plot(bigH'-60)
    % title(obj.musc_obj{musc}.muscle_name(4:end),'Interpreter','none','FontSize',18)
    % ylabel('Input Synergy Voltage (mV)')
    % subplot(4,2,2)
    % plot(syns*1000)
    % ylabel('Synergy Neuron Voltage (mV)')
    % subplot(2,1,1)
    set(gcf,'Position',[-(scW-10)/2 130 scW/2 985])
    axes(ha(1))
        plot(xax,1000.*V_musc(muscnum,:)','LineWidth',2)
        title(obj.musc_obj{muscnum}.muscle_name(4:end),'Interpreter','none','FontSize',18)
        ylabel('MN Activation (mV)')
        hold on
        plot(xax,1000.*keymns(:,mnInd),'LineWidth',2)
        legend({'Calculated';'Simulated'},'Location','northeast')
    %subplot(2,1,2)
    axes(ha(2))
        plot(xax,forces(:,muscnum),'LineWidth',2)
        ylabel('Muscle Tension (N)')
        hold on
        plot(xax,keytens(:,tnInd),'LineWidth',2)
        legend({'Calculated';'Simulated'},'Location','northeast')
    axes(ha(3))
        plot(xax,abs(forces(:,muscnum)-keytens(:,tnInd)),'LineWidth',2)
        ylabel('Muscle Tension (N)')
        legend({'Difference'},'Location','northeast')
        
    %% Passive
    if 0
        [b,ks,kp,Lw,Lr,xoff,Fmax,steepness,mL,mV,ST_max,yoff,dt] = getMuscParams(obj,muscnum,forces);
        T = forces(:,muscnum);
        Tdot = gradient(T,dt);
        delL_musc = max(mL-Lr,0);
                    active = (b/ks).*Tdot+(1+kp/ks).*T;
                passive = kp.*delL_musc+b.*mV;
               figure
                subplot(2,1,1)
                plot(active)
                ylabel('active')
                subplot(2,1,2)
                plot(passive)
                ylabel('passive')
    end
    %% Function: getMuscParams    
    function [b,ks,kp,Lw,Lr,xoff,Fmax,steepness,mL,mV,ST_max,yoff,dt2] = getMuscParams(obj,mnum,forces)
        n2 = 500;
        [beg,ennd] = obj.find_step_indices;
        m2 = length(obj.theta_motion(beg:ennd));
       musc = obj.musc_obj{mnum};
       b = musc.damping;
       ks = musc.Kse;
       kp = musc.Kpe;
       if mnum == 8
        Lw = musc.l_width*1.5;
       else
        Lw = musc.l_width;
       end
       Lr = musc.RestingLength;
      xoff = musc.x_off;
       Fmax = musc.max_force;
       steepness = musc.steepness;
       mL = interp1(1:m2,musc.muscle_length_profile(beg:ennd),linspace(1,m2,n2))';
       mV = interp1(1:m2,musc.muscle_velocity_profile(beg:ennd),linspace(1,m2,n2))';
       ST_max = musc.ST_max;
       yoff = musc.y_off;
       dt2 = ((ennd-beg)*obj.dt_motion)/length(forces);
    end
st_curve = @(Fmax,steepness,xoff,V,yoff) (Fmax./(1+exp(steepness*(xoff-V))))+yoff;
fl = @(Lm,Lr,Lw) 1-(Lm-Lr).^2./Lw^2;
end