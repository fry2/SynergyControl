function F = equilsolver_eqn(x,Lr,Fo)
% This solver is used with equilsolver_func.m to calculate Ks, Kp, and Am values that satisfy a pre-defined FL curve
% Refers to notes taken on 5-13-2019
% Revisions made 2-28-2020

%     x(1) = ks;
%     x(2) = kp;
%     x(3) = Am;   

% % @ L=Lr, F=Fo
%     F(1) = (x(3)./(1+x(2)./x(1)))-Fo;
% % @ L=1.4Lr, F=.36Fo
%     F(2) = ((.36.*x(3)+.4.*x(2).*Lr)./(1+x(2)./x(1)))-.5*Fo;
% % @ L=1.49Lr, F=.0396*Fo
%     F(3) = ((.0396.*x(3)+.49.*x(2).*Lr)./(1+x(2)./x(1)));

% @ L=Lr, F=Fo
   % F(1) = ((1.05*Fo)./(1+x(2)./x(1)))-Fo;
% @ L=1.4Lr, F=.36Fo
    F(1) = ((.36.*(1.05*Fo)+.4.*x(2).*Lr)./(1+x(2)./x(1)))-.5*Fo;
% @ L=1.49Lr, F=.0396*Fo
    F(2) = ((.0396.*(1.05*Fo)+.49.*x(2).*Lr)./(1+x(2)./x(1)));
end