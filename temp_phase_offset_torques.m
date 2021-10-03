% The purpose of this script is to calculate the phase offset of the total joint torque of the hip when the joint is exposed to an oscillating input
% The phase offset is a measure of the impact that viscoelastic and inertial torques have on the total torque
% Inertial torques have an offset of 180, elastic torques have an offset of 0, and viscous torques have an offset of 90
% How do these values change across scales?

inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_HipOsc.asim";
inSimPath = "C:\Users\fry16\OneDrive\Documents\JointDampingOpt\InjectedProject\JointDampingOpt_injected_Standalone.asim";
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_KneeOsc.asim";
lengthVals = exp(linspace(log(.3),log(50),25));
lengthVals = 1;
%lengthVals = 1;
tInfo = zeros(length(lengthVals),4);
joint = find([contains(inSimPath,'Hip'),contains(inSimPath,'Knee'),contains(inSimPath,'Ankle')]);

count = 1;
for ii = 1:length(lengthVals)
    outPath = legScaler(inSimPath,lengthVals(ii)); obj = design_synergy(outPath);
    [total_joint_torques,passive_muscle_torque,inertial_torque] = compute_total_joint_torque(obj,0);
    tnorm = normer(total_joint_torques(:,joint)); jmnorm = normer(obj.theta_motion(:,joint)); 
    pnorm = normer(passive_muscle_torque(:,joint)); inorm = normer(inertial_torque(:,joint));
    [~,jmlocs] = findpeaks(jmnorm); f = 1/(obj.theta_motion_time(jmlocs(3))-obj.theta_motion_time(jmlocs(2)));
    % Total torque
        dt1 = abs(finddelay(jmnorm(285:end),-tnorm(285:end)+1)*.54e-3);
    % VE Torque
        dt2 = abs(finddelay(jmnorm(285:end),-pnorm(285:end)+1)*.54e-3);
    % I Torque
        dt3 = abs(finddelay(jmnorm(285:end),-inorm(285:end)+1)*.54e-3);
    tInfo(ii,:) = [lengthVals(ii),360*f*dt1,360*f*dt2,360*f*dt3];
    disp([num2str(count),' of ',num2str(length(lengthVals))])
    count = count + 1;
end
%%
phi = @(c_0,k_0,m_0,L,L_0,s,T) atand(((c_0.*s.^2.*(L.^3./L_0).*(2.*pi./T))./(k_0.*s.^2.*(L.^3./L_0)+(1/2).*m_0.*9.81.*(L.^4./L_0.^3)-(1/3).*m_0.*(L.^5./L_0.^3).*(2.*pi./T).^2)));

% For finding "sNew"
    obj = design_synergy(inSimPath); moment_output = compute_joint_moment_arms(obj,1,1);
    for ii = 1:38
        mlens(:,ii) = obj.musc_obj{ii}.muscle_length_profile;
    end
    sNew = sum(mean(abs((moment_output(1:38,:)'./1000)./mlens)))/24;
%%
sOrig = sqrt(1e-3); sNew = .2497; s2use = 'new';

switch s2use
    case 'orig'
        sIns = sOrig;
    case 'new'
        sIns = sNew;
end
rat = phi(1.3e3,12.e3,12,.2*lengthVals,1,sIns,1/f); rat(rat<0) = rat(rat<0)+180;

figure('Position',[681,295,1083,684]);
semilogx(.2*tInfo(:,1),tInfo(:,2),'LineWidth',2);hold on;semilogx(.2*lengthVals,rat,'LineWidth',2)
title('Comparing the Phase Offset of the Rat Model, s = .2487','FontSize',16)
xlabel('Length Scale (m)','FontSize',14); ylabel('Phase Offset (deg)','FontSize',14);legend({'Model';'Sutton-Szczecinski'},'FontSize',16,'Location','southeast')
%% Relating B and Kp to Phase Offset
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_HipOsc.asim";
lengthVals = 1;
tInfo2 = zeros(length(lengthVals),4);
ksval = 100;
% 
% bvals = meshgrid(linspace(1,40,10));
% kpvals = zeros(size(bvals));
% for ii = 1:size(bvals,2)
%     maxkp = bvals(1,ii)/.54e-3-ksval;
%     kpvals(:,ii) = exp(linspace(log(.1*bvals(1,ii)),log(maxkp),10));
% end

bvals = logspace(0,1.4771,25);
kpvals = logspace(-3,3.6532,25);
for ii = 1:length(bvals)
    if (bvals(ii)-.54e-3*ksval)/.54e-3 < kpvals(ii)
        kpvals(ii) = (1-eps)*(bvals(ii)-.54e-3*ksval)/.54e-3;
    end
end
phiout = zeros(size(kpvals));

count = 1;
for ii = 1:size(kpvals,1)
    for jj = 1:size(kpvals,2)
        outPath = legScaler(inSimPath,1,[ksval,kpvals(ii,jj),bvals(ii,jj)]); obj = design_synergy(outPath);
        [total_joint_torques,passive_muscle_torque,inertial_torque] = compute_total_joint_torque(obj,0);
        tnorm = normer(total_joint_torques(:,1)); jmnorm = normer(obj.theta_motion(:,1)); 
        pnorm = normer(passive_muscle_torque(:,1)); inorm = normer(inertial_torque(:,1));
        [~,jmlocs] = findpeaks(jmnorm); f = 1/(obj.theta_motion_time(jmlocs(3))-obj.theta_motion_time(jmlocs(2)));
            % Total torque
                dt1 = abs(finddelay(jmnorm(285:end),-tnorm(285:end)+1)*.54e-3);
            % VE Torque
                dt2 = abs(finddelay(jmnorm(285:end),-pnorm(285:end)+1)*.54e-3);
            % I Torque
                dt3 = abs(finddelay(jmnorm(285:end),-inorm(285:end)+1)*.54e-3);
        tInfo2(count,:) = [1,360*f*dt1,360*f*dt2,360*f*dt3];
        if 360*f*dt2 > 180
            phiout(ii,jj) = 360*f*dt2 - 180;
        else
            phiout(ii,jj) = 360*f*dt2;
        end
        phiout(ii,jj) = 360*f*dt2;
        disp([num2str(count),' of ',num2str(size(kpvals,1)*size(kpvals,2))])
        count = count + 1;
    end
end
%%
cm = autumn(180);

figure;surf(kpvals,bvals,phiout);pbaspect([1 1 1]);
set(gca,'xscale','log','yscale','log');xlabel('Kp','FontSize',16);ylabel('B','FontSize',16);zlabel('Phase Offset (deg)','FontSize',16);colormap(gca,cm(1:max(phiout),:))
%% Relating B and Ks to Phase Offset
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_HipOsc.asim";
lengthVals = 1;
tInfo2 = zeros(length(lengthVals),4);

bvals = meshgrid(5*logspace(0,2,10));
ksvals = zeros(size(bvals));
for ii = 1:size(bvals,2)
    maxks = bvals(1,ii)/.54e-3-100;
%     underB = linspace(.1*bvals(1,ii),bvals(1,ii),5);
%     overB = linspace(bvals(1,ii),maxks,5);
    %ksvals(:,ii) = [linspace(.1*bvals(1,ii),.8*bvals(1,ii),5),linspace(bvals(1,ii),maxks,5)];
    ksvals(:,ii) = exp(linspace(log(.1*bvals(1,ii)),log(maxks),10));
    %ksvals(:,ii) = linspace(maxks*.01,maxks,size(bvals,2));
end
phiout = zeros(size(kpvals));

count = 1;
for ii = 1:size(ksvals,1)
    for jj = 1:size(ksvals,2)
        outPath = legScaler(inSimPath,1,[ksvals(ii,jj),100,bvals(ii,jj)]); obj = design_synergy(outPath);
        [total_joint_torques,passive_muscle_torque,inertial_torque] = compute_total_joint_torque(obj,0);
        tnorm = normer(total_joint_torques(:,1)); jmnorm = normer(obj.theta_motion(:,1)); 
        pnorm = normer(passive_muscle_torque(:,1)); inorm = normer(inertial_torque(:,1));
        [~,jmlocs] = findpeaks(jmnorm); f = 1/(obj.theta_motion_time(jmlocs(3))-obj.theta_motion_time(jmlocs(2)));
            % Total torque
                dt1 = abs(finddelay(jmnorm(285:end),-tnorm(285:end)+1)*.54e-3);
            % VE Torque
                dt2 = abs(finddelay(jmnorm(285:end),-pnorm(285:end)+1)*.54e-3);
            % I Torque
                dt3 = abs(finddelay(jmnorm(285:end),-inorm(285:end)+1)*.54e-3);
            tInfo2(ii,:) = [1,360*f*dt1,360*f*dt2,360*f*dt3];
            if 360*f*dt2 > 180
                phiout(ii,jj) = 360*f*dt2 - 180;
            else
                phiout(ii,jj) = 360*f*dt2;
            end
        phiout(ii,jj) = 360*f*dt2;
        disp([num2str(count),' of ',num2str(size(ksvals,1)^2)])
        count = count + 1;
    end
end
%% Relating Kp and Ks to Phase Offset
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_HipOsc.asim";
lengthVals = 1;
tInfo2 = zeros(length(lengthVals),4);

% [kpvals,ksvals] = meshgrid(linspace(.001,1800,5));
[kpvals,ksvals] = meshgrid(1.85.*logspace(-3,3,10));
phiout = zeros(size(kpvals));

count = 1;
for ii = 1:size(ksvals,1)
    for jj = 1:size(ksvals,2)
        if 1/(kpvals(ii,jj)+ksvals(ii,jj)) >= .54e-3
            outPath = legScaler(inSimPath,1,[ksvals(ii,jj),kpvals(ii,jj),1]); obj = design_synergy(outPath);
            [total_joint_torques,passive_muscle_torque,inertial_torque] = compute_total_joint_torque(obj,0);
            tnorm = normer(total_joint_torques(:,1)); jmnorm = normer(obj.theta_motion(:,1)); 
            pnorm = normer(passive_muscle_torque(:,1)); inorm = normer(inertial_torque(:,1));
            [~,jmlocs] = findpeaks(jmnorm); f = 1/(obj.theta_motion_time(jmlocs(3))-obj.theta_motion_time(jmlocs(2)));
            % Total torque
                dt1 = abs(finddelay(jmnorm(285:end),-tnorm(285:end)+1)*.54e-3);
            % VE Torque
                dt2 = abs(finddelay(jmnorm(285:end),-pnorm(285:end)+1)*.54e-3);
            % I Torque
                dt3 = abs(finddelay(jmnorm(285:end),-inorm(285:end)+1)*.54e-3);
            tInfo2(ii,:) = [1,360*f*dt1,360*f*dt2,360*f*dt3];
            if 360*f*dt2 > 180
                phiout(ii,jj) = 360*f*dt2 - 180;
            else
                phiout(ii,jj) = 360*f*dt2;
            end
        else
            phiout(ii,jj) = NaN;
        end
        disp([num2str(count),' of ',num2str(size(ksvals,1)^2)])
        count = count + 1;
    end
end
%%
function outVec = normer(inVec)
    outVec = (inVec - min(inVec)) / (max(inVec) - min(inVec));
end