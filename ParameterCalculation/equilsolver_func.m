function [params,soln] = equilsolver_func(Lr,Fo)
% For a given resting length and optimal muscle force, calculate Kse, Kpe, and Am based on a pre-defined LT relationship (stored in equilsolver_eqn)  
% Input: Lr: optimal resting length in meters
% Input: Fo: optimal force in Newtons
% Output: params: an 1x3 array containing [Kse, Kpe, Am]
% Output: soln: (~)x3 solution array containing other possible solutions to the optimization equation

    % Test data:
        % BFP data in m and N
        %     Lr = .0376;
        %     Fo = 12.49;
    
    % Optimization cost function. Defines three equations for an LT curve relationship
    equilsolver = @(x) equilsolver_eqn(x,Lr,Fo);
    
    rng default % for reproducibility
    N = 100; % try 100 random start points
    pts = 1000*rand(N,2);
    soln = zeros(N,2); % allocate solution
    fval = zeros(N,2); % allocate solution
    opts = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','MaxFunctionEvaluations',1000);
    for k = 1:N
        [soln(k,:),fval(k,:)] = fsolve(equilsolver,pts(k,:),opts); % find solutions
    end
    
    soln = soln(soln(:,2)>0,:);
    bb = sortrows(soln(sum(soln>0,2)==2,:),1,'descend');
%     kpemax = sortrows(soln,2,'descend');
%     params = kpemax(1,:);
    
    %This reduces the complex solution environment to a single result
    %Vary this line to change the solution
    params = soln((soln(:,3) == min(soln(:,3))),:);
    
    %% Calculate a parametric equation for this muscle
%     xyz=soln;
%     r0=mean(xyz);
%     [~,~,V]=svd(bsxfun(@minus,xyz,r0),0);
%     Ksfit=r0(1)+t*V(1,1);
%     Kpfit=r0(2)+t*V(2,1);
%     Amfit=r0(3)+t*V(3,1);
%     t=0;
%     paramvals = r0'+t*V(:,1);
    
    %% Plot the LT curves defined by the full set of solution parameters (usually all very similar)
%    equil_ltcurve = @(x,L,Lr) (1./(1+x(2)./x(1)))*(x(2).*max(0,L-Lr)+x(3)*(1-((L-Lr).^2)./(.5*Lr)^2));
%     
%     L = linspace(.5*Lr,1.5*Lr,1000);
%     
%     for ii = 1:length(soln)
%         bb = equil_ltcurve(soln(ii,:),L,Lr);
%         plot(L,bb)
%         hold on
%     end
end