function [nsys] = indivProjectBuilder(projPath,simPath,forces,obj)

    simPath = char(simPath); projPath = char(projPath);
    muscle_out = scrapeFileForMuscleInfo(simPath);
    numMuscles = length(muscle_out);
    docDir = fileparts(simPath);
    %disp(['Starting to build Animatlab ',docType,' file "',docName,docType,'".'])
    fclose('all');
    delete([docDir,'\Trace*'])
    delete([pwd,'\Data\InjectedCurrent\*'])

    nsys = CanvasModel;

    %% The Time Debacle
    %[beg,ennd] = obj.find_step_indices;
    times = obj.theta_motion_time([obj.sampling_vector(1),obj.sampling_vector(end)]);
    simTime = diff(times);
    %simTime = max(obj.theta_motion_time);
    nsys.proj_params.simendtime = simTime;
    nsys.proj_params.physicstimestep = 1000*(obj.dt_motion); % dt in ms

    [~,V_musc] = Am_generator(obj,forces');
    current2inject = 1000.*(V_musc+.06);
    ci2 = movmean(current2inject',floor(length(current2inject)/3));
%     constvals = mean(ci2);
    
   % current2inject(:,1:3705) = 0;

    %Build a line of motor neurons and muscles first
    for ii = 1:numMuscles
        neurpos = [(ii)*27 250];
        nsys.addItem('n', neurpos, [1000 1000])
        nsys.addMuscle([muscle_out{ii,2},'-neural'],muscle_out{ii,3},neurpos+[-25.9 150])
        nsys.addLink(nsys.neuron_objects(ii),nsys.muscle_objects(ii),'adapter')
        nsys.addStimulus(nsys.neuron_objects(ii),'dc')
        nsys.stimulus_objects(ii).starttime = 0;
        nsys.stimulus_objects(ii).endtime = simTime;
        [outData,stimPath] = generate_direct_current_file(nsys,current2inject(ii,:),nsys.stimulus_objects(ii).name);
        nsys.stimulus_objects(ii).current_wave = outData;
        nsys.stimulus_objects(ii).current_data_file = stimPath;
    end

    % Build datatool viewers for high action muscles to provide insight into motorneuron and muscle activity
    %[muscForces,muscInds] = sortrows(max(forces)',1,'descend');
    %muscles2check = muscInds(1:5);
    muscles2check = 1:length(obj.musc_obj);
    numDTs = size(nsys.datatool_objects,1);
    nsys.addDatatool({'KeyMNs','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscleLen','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscTen','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscAct','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscAl','endtime',nsys.proj_params.simendtime-.01})
    % For a given muscle, find that muscle's adapter information and then find the *neuron* that's feeding that adapter
    for ii = 1:length(muscles2check)
        adInd = find(contains({nsys.adapter_objects(:).destination_node_ID},nsys.muscle_objects(muscles2check(ii)).ID));
        nInd = find(contains({nsys.neuron_objects(:).ID},nsys.adapter_objects(adInd).origin_node_ID));
        nsys.addDTaxes('KeyMNs',nsys.neuron_objects(nInd),'MembraneVoltage')
        nsys.addDTaxes('KeyMuscleLen',nsys.muscle_objects(muscles2check(ii)),'MuscleLength')
        nsys.addDTaxes('KeyMuscTen',nsys.muscle_objects(muscles2check(ii)),'Tension')
        nsys.addDTaxes('KeyMuscAct',nsys.muscle_objects(muscles2check(ii)),'Activation')
        nsys.addDTaxes('KeyMuscAl',nsys.muscle_objects(muscles2check(ii)),'Tl')
    end

    nsys.create_animatlab_simulation(simPath);
    nsys.create_animatlab_project(projPath);
    %disp(['Animatlab ',docType,' file "',docName,docType,' created.'])
%%
    function equation = generate_synergy_eq(bigH,simTime,dt,projBool)
        %dt = .00054;
        time = (0:dt:simTime)';

        if sum(size(bigH)) == 2
            coeffs = cell(3,2);
            coeffs(1:3,1) = [{'a1'};{'b1'};{'c1'}];
            coeffs(1:3,2) = [{bigH}, {0}, {1.5708}];
            equation = sum_of_sines_maker(coeffs,projBool);
        else

            bigH = interpolate_for_time(time,bigH);

            % Create a sum of sines equation for the joint angle waveforms
            fitresult = sumsinesFit(time, bigH,8);
            % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
            coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
            % Equations are in the format necessary for integration into Animatlab's .asim filetype
            equation = sum_of_sines_maker(coeffs,projBool);
        end
    end
%%
    function waveformsBig = interpolate_for_time(time,waveforms)
        % Interpolate the undersampled input to match the required time vector
        waveforms = waveforms';
        m = length(time);
        n = length(waveforms);
        if m ~= n
            waveformsBig = interp1(1:n,waveforms,linspace(1,n,m));
        else
            waveformsBig = waveforms';
        end

        avgblocks = floor(.01*length(waveformsBig));
        coeffblocks = ones(1,avgblocks)/avgblocks;
        for i=1:2
            waveformsBig = filtfilt(coeffblocks,1,waveformsBig);
        end
        waveformsBig = waveformsBig';
    end
end