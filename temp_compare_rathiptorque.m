function temp_compare_rathiptorque()
    obj = design_synergy("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_Size_Rat_Walk_Rat.asim");
    [t_t,t_m,t_i] = compute_total_joint_torque(obj,0);
    results_cell = pedotti_optimization(obj);
    force_mn = results_cell{3,2};
    rm = find_relevant_muscles(obj);
    for kk = 1:3
        [moment_arms] = compute_joint_moment_arms(obj,kk,1);
        moment_arms = moment_arms(rm{kk},:)'./1000;
        active_forces = force_mn(:,rm{kk});
        active_torques(:,kk) = sum(moment_arms.*active_forces,2);
    end
    figure;
    subplot(2,1,1)
    plot(t_t)
    subplot(2,1,2)
    plot(active_torques)
end