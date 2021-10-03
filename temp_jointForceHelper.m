curVars = whos;
if ~any(contains({curVars.name},'obj'))
    obj = design_synergy("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim");
end

jointNum = 3;
jointAngle = -60;

[moment_output] = compute_joint_moment_arms(obj,jointNum,1);

[jointMin,jointMax] = bounds(moment_output(end,:));

while jointAngle < jointMin || jointAngle > jointMax
    jointAngle = input(['Joint Angle ',num2str(jointAngle),' is outside the joint bounds [',num2str(jointMin),' ',num2str(jointMax),']. Please input a valid joint angle: \n']);
end

[~,angleInd] = min(abs(moment_output(end,:)-jointAngle));
numMusc = size(obj.musc_obj,1);
muscBool = sum(moment_output(1:numMusc,:),2)'~=0;
momentArms = moment_output(muscBool,angleInd);

posArms = momentArms>0;
negArms = momentArms<0;


% count = 1;
% lb = zeros(1,length(momentArms));
% ub = lb;
% for ii=1:length(muscBool)
%     if muscBool(ii)
%         lb(count) = 0;
%         ub(count) = obj.musc_obj{ii}.max_force;
%         Fmax(count) = ub(count);
%         count = count + 1;
%     end
% end
% beq = 0;
% x0 = ub;
% fun = @(x) sum((x./Fmax).^2);
% options = optimoptions('fmincon','Display','none','Algorithm','sqp','OptimalityTolerance',1e-4);
% 
% [force,~,exitflag] = fmincon(fun,x0,[],[],momentArms',beq,lb,ub,[],options);

taForce = obj.musc_obj{5}.max_force/2;
taTorque = taForce*moment_output(5,angleInd);
mgfdlRatio = .5;

syms mg fdl

[S] = solve(.5*mg == fdl,mg*moment_output(4,angleInd)+fdl*moment_output(6,angleInd) == -taTorque);

mgForce = double(S.mg);
fdlForce = double(S.fdl);

clear mg fdl

mg = obj.musc_obj{4};
ta = obj.musc_obj{5};
fdl = obj.musc_obj{6};

forces = [mgForce,taForce,fdlForce];

%Will have to find a way to get these automatically. For now I just ran the simulation and copy/pasted the PassiveTension data
passiveForces = [0.1789034 0 0.05849105];

actForces = forces - passiveForces;

syms V

for ii = 4:6
    musc = obj.musc_obj{ii};
    soln = solve(1+exp(musc.steepness*(musc.x_off-V)) == musc.ST_max/actForces(ii-3));
    Iapp(ii-3) = (double(soln)+.06)*1000;
end

fprintf('\n')
for ii=1:length(muscBool)
    if muscBool(ii)
        fprintf([obj.musc_obj{ii}.muscle_name(4:end),' %.2f mm about the %s at %.2f degrees\n']...
            ,moment_output(ii,angleInd),obj.joint_obj{jointNum}.name(4:end),moment_output(end,angleInd))
    end
end