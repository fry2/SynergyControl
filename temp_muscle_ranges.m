sim_file = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim";
[obj,sim_file,joints,bodies,joint_limits,joint_profile] = design_synergy(sim_file);

for ii = 1:length(obj.musc_obj)
    muscle = obj.musc_obj{ii};
    mL = muscle.muscle_length_profile;
    muscle_ranges{ii,1} = muscle.muscle_name;
    if min(mL) < muscle_ranges{ii,2}
        muscle_ranges{ii,2} = min(mL);
    end
    if max(mL) > muscle_ranges{ii,3}
        muscle_ranges{ii,3} = max(mL);
    end
end