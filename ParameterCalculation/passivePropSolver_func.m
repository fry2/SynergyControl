function F = passivePropSolver_func(x,Lr,Fo)
    % Solving for ks, kp, and STmax such that when the muscle is maximally stimulated at steady state, the force produced is less than or equal to Fmax
    % Refers to work in notes from blue 3/6/2020
        F(1) = Fo - (x(1)/(x(1)+x(2)))*(x(3)+(x(2)*Lr)^2/(16*x(3)));
        %For an LT curve with Lr = Lmin, Lwidth = Lmax-Lmin, and beta = 1/(Lmax/Lmin-1)
        % must add a beta parameter to input variables
        %F(1) = Fo - (x(1)/(x(1)+x(2)))*(x(3)+(x(2)^2*Lr^2)/(4*beta^2*x(3)));
end