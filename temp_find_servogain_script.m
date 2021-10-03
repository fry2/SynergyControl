% The purpose of this script is to test how the servo gain parameters (Max Motor Torque, Max Velocity, and Servo Gain) impact motion/passive tension
% at small scalesoutPath = legScaler(inSimPath,.025); inText = importdata(outPath);
inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_walking.asim";
obj = design_synergy(inSimPath);

% Define what a "good" joint angle profile looks like
goodTheta = obj.theta_motion;

factor = .025;
outPath = legScaler(inSimPath,factor); inText = importdata(outPath);

% Max Motor Torque, Max Joint Velocity, Servo Gain
servoGainParams = [10000,100,1000];
[outVal,outObj] = temp_find_servogain(servoGainParams,inText,goodTheta);

%%
figure;subplot(2,1,1);plot(outObj.theta_motion);xlim([0 3055]);title({['Scale = ',num2str(factor),'x Rat Size'];'Joint Motion'},'FontSize',14);
subplot(2,1,2);plot(outObj.passive_tension);xlim([0 3055]);title('Passive Muscle Tension','FontSize',14)