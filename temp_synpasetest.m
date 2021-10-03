syntest = CanvasModel();
proj_file = [pwd,'\Animatlab\SynergyWalking\SynergyWalking20200109.aproj'];
meshMatch(proj_file)
revised_file = strcat(proj_file(1:end-6),'_fake.aproj');

syntest.addItem('n',[125 50],[1000 1000])
syntest.addStimulus(syntest.neuron_objects(1))
t = 0:.001:10;
inEq = 10*sin(t)+20;
    fitresult = sumsinesFit(t,inEq,1);
    % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
    coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
    % Equations are in the format necessary for integration into Animatlab's .asim filetype
    equation = sum_of_sines_maker(coeffs,1);
% syntest.stimulus_objects(1).eq = equation;
syntest.stimulus_objects(1).eq = equations{1};
syntest.stimulus_objects(1).starttime = 0;
syntest.stimulus_objects(1).endtime = 10;

syntest.addDatatool('SynapseTest')
syntest.addDTaxes(syntest.datatool_objects(1),syntest.neuron_objects(1),'MembraneVoltage')
kmat = [1,.5,.25,0];
for ii = 1:4
    syntest.addItem('n',[ii*50 100],[1000 1000])
%     syntest.createSynapseType({['syntest-',num2str(ii+1)],'delE',194,'k',kmat(ii)})
    gsyn = (20*kmat(ii))/(2*194);
    syntest.createSynapseType({['syntest-',num2str(ii+1)],'delE',194,'max_syn_cond',gsyn})
    syntest.addLink(syntest.neuron_objects(1),syntest.neuron_objects(ii+1),syntest.synapse_types(ii+4).name)
    syntest.addDTaxes(syntest.datatool_objects(1),syntest.neuron_objects(ii+1),'MembraneVoltage')
end

syntest.addItem('n',[75 50],[1000 1000])
syntest.addStimulus(syntest.neuron_objects(6))
syntest.stimulus_objects(2).starttime = 0;
syntest.stimulus_objects(2).endtime = 10;
syntest.stimulus_objects(2).eq = equations{2};
gsyn2 = (20*.5)/(2*194);
numSyns = length(syntest.synapse_types);
syntest.createSynapseType({['syntest-',num2str(ii+2)],'delE',194,'max_syn_cond',gsyn2})
syntest.addLink(syntest.neuron_objects(6),syntest.neuron_objects(2),syntest.synapse_types(numSyns+1).name)
syntest.addDTaxes(syntest.datatool_objects(1),syntest.neuron_objects(6),'MembraneVoltage')

%     syntest.synapse_types(5).presyn_thresh = -120;
%     syntest.synapse_types(5).presyn_sat = 0;
%     syntest.synapse_types(numSyns+1).presyn_thresh = -120;
%     syntest.synapse_types(numSyns+1).presyn_sat = -0;

syntest.create_animatlab_project(proj_file);
disp('Finished')