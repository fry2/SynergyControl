function [outParams,jointMotionDeg,outParams_approx,moment_output,mLen,obj,zeta] = muscle_to_joint_VE(inSimPath,joint)
    % A function for calculating joint VE properties from the muscles in the model
    % Input: inSimPath: char or string: path to an ASIM or APROJ file
    % Input: joint: char or double: which joint to focus on (1 - Hip, 2 - Knee, 3, Ankle)
    %% Turn on joint sweeping stimuli and write this to a new sim file
    jointNames = {'Hip';'Knee';'xAnkle'};
    if isnumeric(joint)
        joint = jointNames{joint};
    else
        if ~ischar(joint)
            joint = char(joint);
        end
    end
    
    inText = importdata(inSimPath);
    simEndtime = str2double(extractBetween(string(inText{find(contains(inText,'SimEndTime'))}),'>','</'));
    
    jointId = [contains(joint,'hip','IgnoreCase',1),contains(joint,'knee','IgnoreCase',1),contains(joint,'ankle','IgnoreCase',1)];
    if sum(jointId) ~= 1
        error('Error in joint identity. Joint must identify one: hip, knee, or ankle')
    end
    
    % Disable all joint stimuli just to be safe
    for ii = 1:length(jointNames)
        % Disable MaxMin
            jStimInd = find(contains(inText,['<Name>zMaxMin',jointNames{ii}]));
            enInd = find(contains(inText(jStimInd:end),'<Enabled>'),1,'first')+jStimInd-1;
            inText{enInd} = '<Enabled>False</Enabled>';
        % Disable Constant
            jStimInd = find(contains(inText,['<Name>Constant',jointNames{ii}]));
            enInd = find(contains(inText(jStimInd:end),'<Enabled>'),1,'first')+jStimInd-1;
            inText{enInd} = '<Enabled>False</Enabled>';
        % Disable Walking
            jStimInd = find(contains(inText,['<Name>Walking_',jointNames{ii}]));
            enInd = find(contains(inText(jStimInd:end),'<Enabled>'),1,'first')+jStimInd-1;
            inText{enInd} = '<Enabled>False</Enabled>';
    end
    
    jStimInd = find(contains(inText,['<Name>zMaxMin',jointNames{jointId}]));
    enInd = find(contains(inText(jStimInd:end),'<Enabled>'),1,'first')+jStimInd-1;
    inText{enInd} = '<Enabled>True</Enabled>';
    
    nonJoints = jointNames(~jointId);
    for ii = 1:2
        % Enable Constant
        jStimInd = find(contains(inText,['<Name>Constant',nonJoints{ii}]));
        enInd = find(contains(inText(jStimInd:end),'<Enabled>'),1,'first')+jStimInd-1;
        inText{enInd} = '<Enabled>True</Enabled>';
    end
    
    % Write inText to an ASIM document
    rootPath = fileparts(inSimPath);
    jobSavePath = [char(rootPath),'\muscle_to_joint_VE.asim'];
    %jobSavePath = ['C:\Users\fry16\OneDrive\Documents\temp_files\m2jVE_',jobID,'.asim'];
    fileID = fopen(jobSavePath,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);
    
    obj = design_synergy(jobSavePath);
    
    %% Calculate the joint VE parameters
    theta = obj.joint_obj{jointId}.rec_angle_profile;
    jointShifts = [98.4373 102.226 116.2473];
    %neutralTheta = (90-jointShifts(jointId)).*(pi/180);  
    jointMotionDeg = (180/pi).*theta+jointShifts(jointId);
    %theta = theta - neutralTheta;
    theta = (jointMotionDeg-90).*(pi/180);
    
    mLen = zeros(length(obj.theta_motion),38);
    for ii = 1:length(obj.musc_obj)
        mLen(:,ii) = obj.musc_obj{ii}.muscle_length_profile';
    end
    
    bodyMasses = zeros(1,4); bodyLengths = zeros(1,3);
    for ii = 1:4
        bodyMasses(ii) = obj.body_obj{ii}.mass;
        if ii > 1
            bodyLengths(ii-1) = obj.body_obj{ii}.length;
        end
    end
    
    % Find ROTATIONAL Kgrav
    m2 = 0; 
    if find(jointId) ~= 3
        for ii = find(jointId)+2:4
            m2 = m2 + bodyMasses(ii)/1000; % assumes mass is provided in grams
        end
    end
    nextBod = obj.body_obj{find(jointId)+1}; g = 9.81; L = nextBod.length; m1 = nextBod.mass/1000;
    if find(jointId) == 3
        % have to use toe position for ankle
        bodyVec = (obj.musc_obj{20,1}.pos_attachments{5,4}-obj.joint_obj{jointId}.sim_position_profile)';
    else
        bodyVec = (obj.joint_obj{find(jointId)+1}.sim_position_profile-obj.joint_obj{jointId}.sim_position_profile)';
    end
    vertVec = (obj.joint_obj{jointId}.sim_position_profile+[0,-5,0])';
    if find(jointId) == 2
        % The knee's motion is opposite hip and ankle
        nn = repmat([0;0;-1],1,length(vertVec));
    else
        nn = repmat([0;0;1],1,length(vertVec));
    end
    thetaFromVertical = atan2(dot(nn,cross(vertVec,bodyVec)),dot(vertVec,bodyVec));
    Mgrav = g*L*(m1/2+m2).*sin(thetaFromVertical);
    Kgrav = (Mgrav./thetaFromVertical)';
    Kgrav_approx = repmat(g*L*(m1/2+m2),length(Kgrav),1);
    
    moment_output = compute_joint_moment_arms(obj,find(jointId),1)./1000; % Make sure moment arms are in meters, not mm
    relevant_muscles = moment_output(1:38,1)~=0;
    
    keq = zeros(1,38); beq = keq;
    for ii = 1:38
        musc = obj.musc_obj{ii}; ks = musc.Kse; kp = musc.Kpe;
        keq(ii) = (ks*kp)/(ks+kp);
        beq(ii) = musc.damping;
    end
    
    % Find ROTATIONAL Kelas
    Kelas = (sin(theta)./theta).^2.*sum(keq.*moment_output(1:38,:).^2',2);
    Kelas_approx = sum(keq.*moment_output(1:38,:).^2',2);

    % Find ROTATIONAL Damping, B
    B = sum(beq.*(moment_output(1:38,:)'.*cos(theta)).^2,2);
    B_approx = sum(beq.*(moment_output(1:38,:)').^2,2);
    
    for ii = 1:3
        m(ii) = (1/3).*bodyLengths(ii).^2*sum(bodyMasses(ii+1:4)./1000);
    end
    zeta = B./(2.*sqrt(m(jointId).*(Kelas+Kgrav)));
    
  outParams = [Kelas,B,Kgrav];
  outParams_approx = [Kelas_approx,B_approx,Kgrav_approx];
  
  % Finding positive and negative stiffnesses for a joint, based on moment arms
%   sumArms = sum(moment_output(1:38,:),2);
%   posK = zeros(1,length(moment_output));
%   negK = zeros(1,length(moment_output));
%   posPT = []; negPT = [];
%   for ii = 1:length(sumArms)
%       musc = obj.musc_obj{ii}; lr = musc.RestingLength; pt = musc.passive_tension'; ks = musc.Kse; kp = musc.Kpe;
%       if sumArms(ii) < 0
%           negK = negK + ((ks*kp)/(ks+kp)).*moment_output(ii,:).^2;
%           negPT = [negPT; pt];
%       elseif sumArms(ii) > 0
%           posK = posK + ((ks*kp)/(ks+kp)).*moment_output(ii,:).^2;
%           posPT = [posPT; pt];
%       end
%   end
% figure;
% subplot(4,1,1);plot(obj.theta_motion(:,3).*(180/pi)+116.2473,'LineWidth',3);
% subplot(4,1,2);plot(posK,'LineWidth',3);title(['PosK ',num2str(sum(sumArms>0))]);
% subplot(4,1,3);plot(negK,'LineWidth',3);title(['NegK ',num2str(sum(sumArms<0))])
% subplot(4,1,4);plot(abs(negK)./posK,'LineWidth',3);title('NegK/PosK')
% keyboard
end