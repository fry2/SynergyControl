load('G:\My Drive\Rat\SynergyControl\Data\emgTestSet.mat');

[obj,sim_file,joints,bodies,joint_limits,joint_profile,sdata,passive_tension] = design_synergy(simPath);
%%
% Process data locations 
    docDir = fileparts(simPath);
    fclose('all');
    delete([docDir,'\Trace*'])
    delete([pwd,'\Data\InjectedCurrent\*'])
    delete([pwd,'\Data\ClusterJobs\*'])

% scrape input file for muscle info

    
    scaleVector = [4.2010    0.4372    2.0059    4.7799    4.0554    2.5990    1.6027    1.4557    1.7962    1.2017    3.1493    4.1904];
    lb = zeros(size(scaleVector));
    ub = [4.7399    4.9699    2.1680    2.7364   11.4543    5.3372   25.7270    1.9180    2.9032    3.4796    8.1781   13.4493];
    correctInds = [8     3    12     5     4     9     1    10     2     6     7    11];
    
    for ii = 1:length(ub)
        ub2(ii) = ub(correctInds(ii));
    end
    scaleVector = ub2./2;
    %ub = 5.*ones(size(scaleVector));
    
    %outVal = objFun(scaleVector,emgMns,obj.theta_motion);
    
    options = optimoptions('fmincon','Algorithm','sqp','Display','iter-detailed','MaxFunctionEvaluations',5);
    pattOpts = optimoptions('patternsearch','Display','iter','PlotFcn',{'psplotbestx','psplotbestf','psplotmeshsize'},'MaxTime',2*60,'MaxFunEvals',300,...
        'UseParallel',true);%'MeshContractionFactor',.4,'MeshExpansionFactor',2.5,
    fun = @(scaleVec) objFun(scaleVec,emgMns,obj.theta_motion);
%     [outScale,fVal] = fmincon(fun,scaleVector,[],[],[],[],lb,ub,[],options);
    [outScale,fVal] = patternsearch(fun,scaleVector,[],[],[],[],lb,ub2,[],pattOpts);

    [~,jointMotion] = objFun(outScale,emgMns,obj.theta_motion);
    
figure
plot(jointMotion)
xlabel('Joint Motion')

function [outVal,jointMotion] = objFun(scaleVec,emgCell,normalMotion)
    simPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking_Standalone.asim';
    projPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_EMGwalking.aproj';
    nsys = CanvasModel;
    muscle_out = scrapeFileForMuscleInfo(projPath);
    numMuscles = length(muscle_out);    
    %activation(:,[1:7,9:12]) = 0;
    
    [~,simName] = fileparts(simPath);
    [~,jobID] = fileparts(tempname);
    mkdir([pwd,'\Data\ClusterJobs\',jobID]);
    copyfile(simPath, [pwd,'\Data\ClusterJobs\',jobID]);
    jobSimPath = [pwd,'\Data\ClusterJobs\',jobID,'\',jobID,'.asim'];
    movefile([pwd,'\Data\ClusterJobs\',jobID,'\',simName,'.asim'], [pwd,'\Data\ClusterJobs\',jobID,'\',jobID,'.asim']);
    
    inText = importdata(jobSimPath);
    txtInd = find(contains(inText,'OutputFilename>Joint'));
    inText{txtInd} = replaceBetween(inText{txtInd},'ion','.txt',['_',jobID]);
    fileID = fopen(jobSimPath,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);
    
    activation = emg2activation(emgCell{1});
    
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
            emgColumn = find(contains(emgCell{2}(2,:),muscle_out{ii,2}));
            emgTracker(ii) = emgColumn;
            dataWave = scaleVec(ii)*[activation(:,emgColumn);activation(:,emgColumn);activation(:,emgColumn)]*20;
            dcStimPath = [pwd,'\Data\ClusterJobs\',jobID,'\',jobID,'_',nsys.stimulus_objects(ii).name,'.txt'];
            [outData,stimPath] = generate_direct_current_file(nsys,dataWave,nsys.stimulus_objects(ii).name,dcStimPath);
            nsys.stimulus_objects(ii).current_wave = outData;
            nsys.stimulus_objects(ii).current_data_file = stimPath;
    end

    % build datatools
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

    %     cm = floor(jet(numMuscles)*255);
    %     % recolor axis lines to make it easier to differentiate plot elements
    %     for ii = 1:length(nsys.datatool_objects)
    %         for jj = 1:length(nsys.datatool_objects(ii).axes_objects)
    %             nsys.datatool_objects(ii).axes_objects(jj).linecolor = rgb2anim(cm(jj,:));
    %         end
    %     end

        % write object to proj and sim files

        %savePath = [pwd,'\Data\ClusterJobs\','Job',num2str(w.ProcessId),'\temp.asim'];
        %savePath = [tempname([pwd,'\Data\ClusterJobs\']),'.asim'];
        jobSavePath = [jobSimPath(1:end-5),'_mod',jobSimPath(end-4:end)];
        nsys.create_animatlab_simulation(jobSimPath,jobSavePath);
        %nsys.create_animatlab_project(projPath);

    %%    
    %modSimPath = [simPath(1:end-5),'_fake',simPath(end-4:end)];
    %simTest = processSimData(jobSavePath);
            sour_folder = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin';
            executable = [string([sour_folder,'\AnimatSimulator']),string(jobSavePath)];
            %executable = ['"',sour_folder,'\AnimatSimulator','" "',trialFileName,'"'];
            %[~,~,err] = jsystem(executable,'noshell');
            jsystem(executable);

    %%
    jointMotionPath = [pwd,'\Data\ClusterJobs\',jobID,'\','JointMotion_',jobID,'.txt'];
    ds = importdata(jointMotionPath);
    jointMotion = ds.data(:,3:5);
    %jointMotion = simTest(1).data;
    jointGrades = vecnorm(normalMotion-jointMotion,2,1);
    outVal = 10*jointGrades(1)+10*jointGrades(2)+jointGrades(3);
    
    delete([pwd,'\Data\ClusterJobs\',jobID,'\*'])
    rmdir([pwd,'\Data\ClusterJobs\',jobID])
end