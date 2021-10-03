function C = axisAngleRotation(angle,joint_axis)
    % For na input angle and an input joint axis, find the rotation matrix
    % Input: angle (double): angle in RAD
    % Input: joint_axis (double): 1x3 unit vector of a joint axis
    c = cos(angle);
    s = sin(angle);

    a1 = joint_axis(1);
    a2 = joint_axis(2);
    a3 = joint_axis(3);

    C = [c+a1^2*(1-c), a1*a2*(1-c)-a3*s, a1*a3*(1-c)+a2*s;...
         a1*a2*(1-c)+a3*s, c+a2^2*(1-c), a2*a3*(1-c)-a1*s;...
         a1*a3*(1-c)-a2*s, a2*a3*(1-c)+a1*s, c+a3^2*(1-c)];
end