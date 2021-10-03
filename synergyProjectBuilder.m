projPath = [pwd,'\Animatlab\SynergyWalking\SynergyControl.aproj'];
meshMatch(projPath);
revised_file = strcat(projPath(1:end-6),'_fake.aproj');
[projDir,projName,ext] = fileparts(revised_file);
disp(['Starting to build Animatlab project ',projName,ext])
delete([projDir,'\Trace*'])
%%
original_text = importdata(projPath);
muscle_info = scrapeFileForMuscleInfo(projPath);
numMuscles = length(muscle_info);
numZones = 6;
muscZones = zoning_sorter(original_text,numZones);
muscle_info(:,4) = muscZones(:,2); clear muscZones

nsys = CanvasModel;
neurpos = [];

%% The Time Debacle
simTime = 10;
nsys.proj_params.simendtime = simTime+.01;
nsys.proj_params.physicstimestep = .54; % dt in ms

%Build a line of motor neurons and muscles first
for ii = 1:numMuscles
    neurpos = [(ii)*26 250];
    nsys.addItem('n', neurpos, [1000 1000])
    nsys.addMuscle([muscle_info{ii,2},'-neural'],muscle_info{ii,3},neurpos+[-25.9 150])
    nsys.addLink(nsys.neuron_objects(ii),nsys.muscle_objects(ii),'adapter')
end

%% Create the synergy neurons
% Generate color palette for synergy nodes
% cm = jet(3*numZones);
% cm = round(cm(3:3:end,:)*255);
    cm = round([0,1,0;...
          1,0,0;...
          0.9216,0.8235,0.2039;...
          0,0,1;...
          1,0,1;...
          0,1,1].*255);

syngap = (958.5-numZones*(25))/(numZones+1);
nsys.addDatatool('SynergyStim');

for ii = 1:numZones
    synpos = [(ii-1)*25+(ii)*syngap 100];
    nsys.addItem('n', synpos, [1000 1000]);
    colordec = rgb2anim(cm(ii,:));
    nsys.neuron_objects(numMuscles+ii).color = colordec;
    nsys.neuron_objects(numMuscles+ii).name = ['Syn ',num2str(ii)];
    nsys.neuron_objects(numMuscles+ii).nsize = [100,100];
    nsys.addDTaxes(nsys.datatool_objects(1),nsys.neuron_objects(numMuscles+ii),'MembraneVoltage');
    nsys.addStimulus(nsys.neuron_objects(numMuscles+ii),'tc')
    nsys.stimulus_objects(ii).starttime = 0;
    nsys.stimulus_objects(ii).endtime = simTime;
end

%% Connect synergy neurons to correct MNs
for ii = 1:numMuscles
    mn = mn_from_musc_name(nsys,muscle_info{ii,2});
    zoneNum = muscle_info{ii,4};
    synNeur = nsys.neuron_objects(zoneNum+numMuscles);
            nsys.addLink(synNeur,mn,'SignalTransmission1')
            synapseName = ['Syn-',num2str(zoneNum+numMuscles),'-',num2str(ii)];
            nsys.createSynapseType({synapseName,'delE',194,'k',1})
            numLinks = size(nsys.link_objects,1);
            nsys.link_objects(numLinks).synaptictype = synapseName;
end

% Build datatool viewers for high action muscles to provide insight into motorneuron and muscle activity
%[muscForces,muscInds] = sortrows(max(forces)',1,'descend');
muscles2check = 1:38;%muscInds(1:5);
%muscles2check = 1:38;
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

nsys.create_animatlab_project(projPath);
disp(['Animatlab project file ',projName,ext,' created.'])

function equation = generate_synergy_eq(bigH,simTime)
    dt = .00054;
%     time = (99*dt:dt:10.01-10*dt)';
    time = (0:dt:simTime)';
    
    bigH = interpolate_for_time(time,bigH);
    
    % Create a sum of sines equation for the joint angle waveforms
    fitresult = sumsinesFit(time, bigH,8);
    % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
    coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
    % Equations are in the format necessary for integration into Animatlab's .asim filetype
    equation = sum_of_sines_maker(coeffs,1);
end

function waveformsBig = interpolate_for_time(time,waveforms)
    % Interpolate the undersampled input to match the required time vector
    waveforms = waveforms';
    m = length(time);
    n = length(waveforms);
    if m ~= n
        waveformsBig = interp1(1:n,waveforms,linspace(1,n,m));
    end

    avgblocks = floor(.01*length(waveformsBig));
    coeffblocks = ones(1,avgblocks)/avgblocks;
    for i=1:2
        waveformsBig = filtfilt(coeffblocks,1,waveformsBig);
    end
    waveformsBig = waveformsBig';
end

function mn = mn_from_musc_name(nsys,muscName)
    adInd = find(contains({nsys.adapter_objects.destination_node_ID},muscName));
    mnID = nsys.adapter_objects(adInd).origin_node_ID;
    neurInd = find(contains({nsys.neuron_objects.ID},mnID));
    mn = nsys.neuron_objects(neurInd);
end