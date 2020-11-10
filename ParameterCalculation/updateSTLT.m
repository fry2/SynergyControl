% For a project with correct viscoelastic properties, inject correct ST/LT values

sim_file = [pwd,'\Animatlab\SynergyWalking\SynergyControl_Standalone.asim'];
obj = design_synergy(sim_file);

load([pwd,'\Data\minmaxlrs.mat'],'outData')

inParams = cell(38,1);

for ii = 1:38
    muscle = obj.musc_obj{ii};
    STmax = (muscle.Kse+muscle.Kpe)/muscle.Kse*muscle.max_force;
    % LT info
    LTind = find(contains(outData.data(:,1),muscle.muscle_name(4:end)));
        ltinfo = outData.data(LTind,:);
        Lr = ltinfo{2};
        Lw = ltinfo{3};
        Lmin = ltinfo{4};
        Lmax = ltinfo{5};
        Perc = ltinfo{6}; % Percentage at bounds
    inParams{ii,1} = muscle.muscle_name;
    inParams{ii,2} = muscle.Kse;
    inParams{ii,3} = muscle.Kpe;
    inParams{ii,4} = STmax;
    [steepness,yoff] = stPropSolver(STmax,muscle.max_force);
    inParams{ii,5} = steepness;
    inParams{ii,6} = 1000*yoff;
    [inParams{:,7}] = deal(-60);
    [inParams{:,8}] = deal(-40);
    [inParams{:,9}] = deal(0);
    inParams{ii,10} = STmax;
    inParams{ii,11} = 100*(Lmin); % 'LT<LowerLimit'
    inParams{ii,12} = 100*Lr; % '<RestingLength'
    inParams{ii,13} = 100*(Lmax); % 'LT<UpperLimit'
    inParams{ii,14} = 100*(Lw); % '<Lwidth'
    inParams{ii,15} = muscle.max_force; % '<MaximumTension'
    inParams{ii,16} = Perc; % 'LT<LowerOutput'
    inParams{ii,17} = Perc; % 'LT<UpperOutput'
end

parameter_injector('G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl.aproj',inParams)



