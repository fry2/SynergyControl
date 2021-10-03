function jointMotion = emg_to_walking_simproj(inSimPath,inEMGdata,inScale)
    % Make a simulation file with the input EMG data
    % Input: inSimPath: (char) path to simulation file to copy
    % Input: inEMGdata: (cell) {1} contains EMG data, {2} contains muscle names
    % Input: inScale: (double) vector of length muscle_names that scale the EMG terms

%     inSimPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking_Standalone.asim';
%     inEMGdata{1} = emgData(:,:,2); inEMGdata{2} = emgNames; inScale = outScale(2,:);
%     jointMotion = emg_to_walking_simproj(inSimPath,inEMGdata,inScale);
    
    simPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking_Standalone.asim';
    projPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking.aproj';
    %simPath = 'G:\My Drive\Rat\SynergyControl\EMG_control\SynergyWalking_EMGwalking_Standalone.asim';
    %projPath = 'G:\My Drive\Rat\SynergyControl\EMG_control\SynergyWalking_EMGwalking.aproj';
    nsys = CanvasModel;
    muscle_out = scrapeFileForMuscleInfo(inSimPath);
    numMuscles = length(muscle_out);    
    if nargin < 3
        inScale = ones(1,numMuscles);
    end
    
    inEMGdata = order_emg_for_sim(inSimPath,inEMGdata);
    activation = emg2activation(inEMGdata{1});
    
    for ii=1:numMuscles
        % make neuron
            neurpos = [(ii)*25.90 250];
            nsys.addItem('n', neurpos, [1000 1000]);
        % make muscle
            nsys.addMuscle([muscle_out{ii,2},'-neural'],muscle_out{ii,3},neurpos+[-25.9 150]);
        % make link neuron and muscle
            nsys.addLink(nsys.neuron_objects(ii),nsys.muscle_objects(ii),'adapter')
        % make direct current stim
            nsys.addStimulus(nsys.neuron_objects(ii),'dc')
        % define stimulus waveform
            dataWave = inScale(ii)*activation(:,ii)*20;
            dcStimPath = [pwd,'\Data\InjectedCurrent\',nsys.stimulus_objects(ii).name,'.txt'];
            outData = generate_direct_current_file(nsys,dataWave,dcStimPath);
            nsys.stimulus_objects(ii).current_wave = outData;
            nsys.stimulus_objects(ii).current_data_file = dcStimPath;
    end

    %build datatools
    nsys.addDatatool({'KeyMNs','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscleLen','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscTen','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscAct','endtime',nsys.proj_params.simendtime-.01})
    nsys.addDatatool({'KeyMuscAl','endtime',nsys.proj_params.simendtime-.01})
    % For a given muscle, find that muscle's adapter information and then find the *neuron* that's feeding that adapter
    for ii = 1:numMuscles
        adInd = find(contains({nsys.adapter_objects(:).destination_node_ID},nsys.muscle_objects(ii).ID));
        nInd = find(contains({nsys.neuron_objects(:).ID},nsys.adapter_objects(adInd).origin_node_ID));
        nsys.addDTaxes('KeyMNs',nsys.neuron_objects(nInd),'MembraneVoltage')
        nsys.addDTaxes('KeyMuscleLen',nsys.muscle_objects(ii),'MuscleLength')
        nsys.addDTaxes('KeyMuscTen',nsys.muscle_objects(ii),'Tension')
        nsys.addDTaxes('KeyMuscAct',nsys.muscle_objects(ii),'Activation')
        nsys.addDTaxes('KeyMuscAl',nsys.muscle_objects(ii),'Tl')
    end

        cm = floor(jet(numMuscles)*255);
        % recolor axis lines to make it easier to differentiate plot elements
        for ii = 1:length(nsys.datatool_objects)
            for jj = 1:length(nsys.datatool_objects(ii).axes_objects)
                nsys.datatool_objects(ii).axes_objects(jj).linecolor = rgb2anim(cm(jj,:));
            end
        end

        %write object to proj and sim files
        outSimPath = [inSimPath(1:end-5),'_mod',inSimPath(end-4:end)];
        outProjPath = [outSimPath(1:end-20),'_mod.aproj'];
        nsys.create_animatlab_simulation(inSimPath,outSimPath);
        nsys.create_animatlab_project(projPath);

    %%    
    %modSimPath = [simPath(1:end-5),'_fake',simPath(end-4:end)];
    simTest = processSimData(outSimPath);
    
    jointMotion = simTest(contains({simTest.name},'JointMotion')).data;
end