function [newTheta,newTime] = make_theta_motion_init0(obj)
    % This function's purpose is to take Fischer's original joint motion and (1) shift it so that the hip initializes at 0
    % (2) creates a linear interpolation for the knee and ankle such that they also initialize at zero
    % Initializing joint motion at zero may alleviate some issues from calculating muscle stimuli, since muscles need to "snap" the leg to a different
    % initial condition
    
    originalTheta = obj.theta_motion;
    hipInit0Bnds = 2269:8448;
    threeSteps = [originalTheta(hipInit0Bnds,:);originalTheta(hipInit0Bnds,:);originalTheta(hipInit0Bnds,:);];
    
    pt2meet = 1050; % the point where the knee and ankle waves stop interpolating and begin the actual motion
    jointTheta = @(x,jnum) (threeSteps(pt2meet,jnum)/pt2meet).*x;
    newTheta = [threeSteps(:,1),[jointTheta(0:pt2meet,2)';threeSteps(pt2meet+2:end,2)],[jointTheta(0:pt2meet,3)';threeSteps(pt2meet+2:end,3)]];
    
    newTime = linspace(0,10,length(newTheta))';
    
%     for j = 1:3
%         [fitresult] = sumsinesFit(newTime, newTheta(:,j));
%         % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
%         coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
%         % Equations are in the format necessary for integration into Animatlab's .asim filetype
%         equations{j} = sum_of_sines_maker(coeffs,0);
%         equations_proj{j} = sum_of_sines_maker(coeffs,1);
%     end
end