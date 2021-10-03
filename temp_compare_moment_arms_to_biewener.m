% Compare moment arms across scales to a figure in Biewener 1990, "Biomechanics of Mammalian...", Figure 2

% A walking simulation path
inSimPath = "C:\Users\fry16\OneDrive\Documents\JointDampingOpt\InjectedProject\JointDampingOpt_injected_Standalone.asim";

% Stance Inds, indices for one stance cycle
stanceInds = 1037:1689;
% Biewener data is for "each joint during the middle third of the support phase"
si_lb = stanceInds(1)+floor(.333*length(stanceInds));
si_ub = stanceInds(end)-floor(.333*length(stanceInds));

% Mass and length vectors. Biewener works in body masses, model works on factors of length. We use Suttinski mass scaling rule of m = 12*(L)^3
masses = logspace(-2,3,10);
lengths = (masses./12).^(1/3);
factor = lengths./.2;

% Pre-allocate for speed c
meanArms = zeros(length(factor),3);

count = 1;
for ii = 1:length(factor)
    outPath = legScaler(inSimPath,factor(ii)); obj = design_synergy(outPath);
    for jj = 1:3
        [moment_output] = compute_joint_moment_arms(obj,jj,1);
        for kk = 1:38
            svals(:,kk) = (moment_output(kk,si_lb:si_ub)'./1000)./obj.musc_obj{kk}.muscle_length_profile(si_lb:si_ub);
        end
        % Moment arm data is in (mm), divide by 10 to put it in (cm) like Biewener
        meanArms(ii,jj) = (1/10)*mean(abs(moment_output(sum(moment_output(1:38,:),2)'~=0,si_lb:si_ub)),'all');
        meanS(ii,jj) = mean(abs(svals(:,sum(svals)~=0)),'all');
    end
    musc = obj.musc_obj{28};
    muscInfo{ii,1} = musc.passive_tension; muscInfo{ii,2} = musc.muscle_length_profile; muscInfo{ii,3} = musc.muscle_velocity_profile;
    muscInfo{ii,4} = musc.Kse; muscInfo{ii,5} = musc.Kpe; muscInfo{ii,6} = musc.damping; muscInfo{ii,7} = musc.max_force;
    disp([num2str(count),' out of ',num2str(length(factor))])
    count = count + 1;
end

%% Plot results
figure('Position',[1009,442,912,549]);
loglog(masses,meanArms,'LineWidth',3);hold on;loglog(masses,r,'LineWidth',3);grid on
xlabel('Body Mass (kg)','FontSize',16);ylabel('Moment Arm Length (cm)','FontSize',16)
title('Moment Arm v Body Length','FontSize',18)
legend({'Hip';'Knee';'Ankle';'Biewener'},'Location','southeast','FontSize',16)