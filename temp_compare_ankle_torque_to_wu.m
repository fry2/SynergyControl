function outVal = temp_compare_ankle_torque_to_wu(inVals)   
    % Optimization says that ksval = 0.0141 and kpval = 1.2132 bring model and data closest in line
    %% Optional: Mess with Ks values of the sim file first
    inSimPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_MinMaxAnkle.asim';
    %inSimPath = 'C:\Users\fry16\OneDrive\Documents\JointDampingOpt\InjectedProject\JointDampingOpt_injected_Standalone.asim';
    inText = importdata(inSimPath);
    ksinds = find(contains(inText,'<Kse>')); kpinds = find(contains(inText,'<Kpe>'));
    for ii = 1:length(ksinds)
        % Ks
        currVal = double(extractBetween(string(inText{ksinds(ii)}),'>','</'));
        inText{ksinds(ii)} = replaceBetween(inText{ksinds(ii)},'>','</',num2str(inVals(1)*currVal));
        % Kp
        currVal = double(extractBetween(string(inText{kpinds(ii)}),'>','</'));
        inText{kpinds(ii)} = replaceBetween(inText{kpinds(ii)},'>','</',num2str(inVals(2)*currVal));
    end
        % Write to file
        jobSavePath = [inSimPath(1:end-5),'_mod.asim'];
        fileID = fopen(jobSavePath,'w');
        fprintf(fileID,'%s\n',inText{:});
        fclose(fileID);
    %
    mma = design_synergy(jobSavePath);
    [pjt,passive_joint_motion,passive_muscle_torque] = compute_passive_muscle_torque(mma,0);
    jm = mma.theta_motion.*(180/pi) + [98.4373 102.226 116.2473];
    stepInds = 1581:2145;
    % Data from "Passive elastic properties of the rat ankle" Wu 2012
    wudat = [-36.18,-.24;-18,.22;0,1.16;8.92,1.93;17.75,2.82];
    wudatBig = [interp1(1:length(wudat),wudat(:,1),linspace(1,length(wudat),length(pjt(stepInds,3))))',interp1(1:length(wudat),wudat(:,2),linspace(1,length(wudat),length(pjt(stepInds,3))))'];
    modeldat = [jm(stepInds,3)-90,1000*pjt(stepInds,3)];
    outVal = sum((-wudatBig(:,2)-modeldat(:,2)).^2);
%     %
    figure;
    plot(modeldat(:,1),modeldat(:,2),'LineWidth',3);hold on;
    plot(wudatBig(:,1),-wudatBig(:,2),'LineWidth',3);
    xlim([-40 20]);legend({'Model';'Wu'});title('Passive Joint Torque')
    xlabel('Ankle Joint Angle (deg)'); ylabel('Passive Joint Torque (mN-m)')
% 
%     %% Calculating Animatlabs bounds to get ankle to move from -40 to 20
%     wubnds = (90-[-40 20]-116.2473).*(pi/180);
%     tt = 0:.54e-3:1.7;
%     cc = @(t) mean(wubnds)+diff(wubnds)/2*sin(6*t);
%     figure;plot(tt,cc(tt))
end