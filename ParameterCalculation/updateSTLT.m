% For a project with correct viscoelastic properties, inject correct ST/LT values

sim_file = [pwd,'\Animatlab\SynergyWalking\SynergyControl_Standalone.asim'];
obj = design_synergy(sim_file);

%load([pwd,'\Data\minmaxlrs_old.mat'],'outData')
load([pwd,'\Data\neutral_lengths.mat'],'neutral_lengths')

load([pwd,'\Data\max_forces.mat'],'max_forces')
mfMult = 1.5;

inParams = cell(38,1);

for ii = 1:38
    muscle = obj.musc_obj{ii};
    MFind = find(contains(max_forces.max_forces(:,1),muscle.muscle_name(4:end)));
    STmax = (1+muscle.Kpe/muscle.Kse)*mfMult*max_forces.max_forces{MFind,2}; 
    steepness = 459.5;
    yoff = -.01*STmax;
    % LT info
    LTind = find(contains(neutral_lengths.data(:,1),muscle.muscle_name(4:end)));
        ltinfo = neutral_lengths.data(LTind,:);
        Lr = ltinfo{2};
        Lw = ltinfo{3};
        Lmin = ltinfo{4};
        Lmax = ltinfo{5};
        Perc = ltinfo{6}; % Percentage at bounds
    inParams{ii,1} = muscle.muscle_name;
    inParams{ii,2} = muscle.Kse;
    inParams{ii,3} = muscle.Kpe;
    inParams{ii,4} = STmax;
    %[steepness,yoff] = stPropSolver(STmax,muscle.max_force);
    inParams{ii,5} = steepness;
    inParams{ii,6} = yoff;
    [inParams{:,7}] = deal(-60);
    [inParams{:,8}] = deal(-40);
    [inParams{:,9}] = deal(0);
    inParams{ii,10} = STmax;
    inParams{ii,11} = 100*(Lmin); % 'LT<LowerLimit'
    inParams{ii,12} = 100*Lr; % '<RestingLength'
    inParams{ii,13} = 100*(Lmax); % 'LT<UpperLimit'
    inParams{ii,14} = 100*(Lw); % '<Lwidth'
    inParams{ii,15} = mfMult*max_forces.max_forces{MFind,2}; % '<MaximumTension'
    inParams{ii,16} = Perc; % 'LT<LowerOutput'
    inParams{ii,17} = Perc; % 'LT<UpperOutput'
end

parameter_injector('G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl.aproj',inParams)
