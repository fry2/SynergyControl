simpath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim";
muscleSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\muscleStim.asim";
motorSimPath ="G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\motorStim.asim";

[motorObj,sim_file,joints,bodies,joint_limits,joint_profile] = design_synergy(motorSimPath);
muscleStruct = processSimData(muscleSimPath);
simLengths = muscleStruct(4).data;
clear motormls
for ii = 1:38
    motormls(ii,:) = motorObj.musc_obj{ii}.muscle_length_profile;
end

muscSimNames = muscleStruct(4).data_cols;
oldSimData = simLengths;
newSimData = zeros(size(oldSimData));
for ii = 1:38
    muscObjName = motorObj.musc_obj{ii}.muscle_name;
    simInd = find(contains(muscSimNames,muscObjName(4:end)));
    if isempty(simInd)
        keyboard
    end
    newSimData(:,ii) = oldSimData(:,simInd);
end
simLengths = newSimData;

meanSim = mean(simLengths);
meanCalc = mean(motormls');

figure
a = subplot(3,1,1);
bar(meanSim.*1000);
b = subplot(3,1,2);
bar(meanCalc.*1000);
c = subplot(3,1,3);
bar((meanSim-meanCalc).*1000);

[~,mnum] = max((meanSim-meanCalc).*1000);
figure
subplot(2,1,1)
plot(motormls(mnum,1:end-10))
ylabel('MOTOR')
subplot(2,1,2)
plot(simLengths(1:end-10,mnum))
ylabel('MUSCLE')

[muscForces,muscInds] = sortrows((meanSim-meanCalc)'.*1000,1,'descend');
muscles2check = muscInds(1:5);
for ii = 1:length(muscles2check)
    fprintf([obj.musc_obj{muscles2check(ii)}.muscle_name,'\n'])
end