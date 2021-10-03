function [elasK,gravK,bRot,zeta,elasK_approx,gravK_approx,bRot_approx,obj,jm] = muscle_to_joint_VE_all(inSimPath)
    [temph,hjm,tempha,~,~,~,zeta(:,1)] = muscle_to_joint_VE(inSimPath,1);
    [tempk,kjm,tempka,~,~,~,zeta(:,2)] = muscle_to_joint_VE(inSimPath,2);
    [tempa,ajm,tempaa,~,~,~,zeta(:,3)] = muscle_to_joint_VE(inSimPath,3);
    % Save elastic stiffness
        elasK = [temph(:,1), tempk(:,1), tempa(:,1)];
        elasK_approx = [tempha(:,1), tempka(:,1), tempaa(:,1)];
    % Save damping
        bRot = [temph(:,2), tempk(:,2), tempa(:,2)];
        bRot_approx = [tempha(:,2), tempka(:,2), tempaa(:,2)];
    % Save gravitational stiffness
        gravK = [temph(:,3), tempk(:,3), tempa(:,3)];
        gravK_approx = [tempha(:,3), tempka(:,3), tempaa(:,3)];
        jm = [hjm, kjm, ajm];
        obj = design_synergy(inSimPath);
end