    close all

    file_dir = fileparts(mfilename('fullpath'));
    load([file_dir,'\Data\processedHindlimbAngles.mat'],'BackMean','BackRaw','completeWaves')
    addpath(genpath(pwd))

    trial = 2;
%     waveform = [BackRaw(:,trial,1)-98,BackRaw(:,trial,2)-90,BackRaw(:,trial,3)-116];
    waveform = [completeWaves(:,trial,1)-98,completeWaves(:,trial,2)-90,completeWaves(:,trial,3)-116];
        initAngs = [-0.9656   -0.1941   -0.6678];
        naturalStopping = [-.748 -.01337 -.47064];
        toeoff = [.1583246 -.15 -.562388];
    %waveform = repmat(toeoff,10000,1);
    clear trial
    
    tstart = tic;
    clear obj
    [obj] = jointMotionInjector(motorSim,waveform,0);
    numMuscles = length(obj.musc_obj);
    
    telapsed = toc(tstart);
    if telapsed>60
        mins = num2str(floor(telapsed/60));
        sec = num2str(round(mod(telapsed,60)));
    else
        mins = num2str(0);
        sec = num2str(round(telapsed,2));
    end
    disp(['Waveforms Injected and Simulated.',' (',mins,'m ',sec,'s)'])
    clear mins sec telapsed tstart
    
    % Optimize forces to meet torque demands
    results_cell = obj.pedotti_optimization;
    oforces = results_cell{2,2}';
    if size(oforces,2) ~= numMuscles
        oforces = oforces';
    end

    forces = smoothdata(oforces,'gaussian',20);
    forces(forces<0) = 0;
    pts = zeros(numMuscles,length(obj.theta_motion));
    
    for ii = 1:numMuscles
        pts(ii,:) = obj.musc_obj{ii}.passive_tension;
    end
    
    %Subtract out the passive portion of the forces
    forces = oforces-pts(:,obj.sampling_vector)';
    
    synergysort = 'currents';
    switch synergysort
        case 'forces'
            [r2scores,recompiled,W,H] = NMFdecomposition(5,forces,0,.04);

            % Coefficients of similarity between original forces and recompiled forces from synergies
            coeffMat = pearsonTest(forces,recompiled,0);
    
            % Am signals that will generate the desired forces
            [Am_musc,V_musc] = Am_generator(obj,forces');
            current2inject = 1000.*(V_musc+.06);
            wave2plot = forces;
        case 'currents'
            % Am signals that will generate the desired forces
            [Am_musc,V_musc,Al_musc_all] = Am_generator(obj,forces');
            current2inject = 1000.*(V_musc+.06);
%             [r2scores,recompiled,W,H] = NMFdecomposition(5,current2inject',0,.04);
            [r2scores,recompiled,W,H] = NMFdecomposition(6,current2inject',0,0);
            wave2plot = current2inject';   
    end
    
    if isfile([file_dir,'\Data\h_equations.mat'])
        delete([file_dir,'\Data\h_equations.mat'])
    end
    
    if isfile([file_dir,'\Data\indiv_equations.mat'])
        delete([file_dir,'\Data\indiv_equations.mat']) 
    end
    
%     figure('name','InjectedWaveforms')
%     plot(obj.theta_motion)
%     title('Joint Angle Waveforms')
%     legend({obj.joint_obj{1}.name(4:end),obj.joint_obj{2}.name(4:end),obj.joint_obj{3}.name(4:end)},'Interpreter','none')
%     
%     figure('name','Current2Inject')
%     plot(current2inject')
%     title('Inject this current into sim')
%     ylabel('Current (nA)')
%     
%     figure('name','OptimizedVSynergies')
%     subplot(2,1,1)
%     plot(wave2plot)
%     title('Optimized Necessary Forces v. Forces Recompiled from Synergies')
%     xlabel('Optimized Muscle Forces')
%     ylabel('Forces(N)')
%     subplot(2,1,2)
%     plot(recompiled)
%     xlabel('Forces Recompiled from Synergies')
%     ylabel('Forces(N)')
    
    %plotWH(wave2plot,W,H,0)