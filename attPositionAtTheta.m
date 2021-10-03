function outPos = attPositionAtTheta(obj,theta,muscAtt)
    % This function is used to determine the global position of a point on the hindlimb for any input leg configuration, theta
    % Input: theta (double): 1x3 vector of joint angles in RAD
    %Input: muscAtt: 1x4 cell array with information about muscle attachment point
        % [attachment intial position, attachment name, body number of attachment, attachment profile for input waveform]
    for ii = 1:3
        limBool = [theta(ii) > max(obj.joint_obj{ii}.limits) theta(ii) < min(obj.joint_obj{ii}.limits)];
        if any(limBool) && obj.joint_obj{ii}.enable_limit
            warning(['Desired theta ',num2str(theta(ii)),' for ',obj.joint_obj{ii}.name,' is outside joint limits, '...
                '[',num2str(obj.joint_obj{ii}.limits(1)),', ',num2str(obj.joint_obj{ii}.limits(2)),']'])
            theta(ii) = obj.joint_obj{ii}.limits(limBool);
        end
    end

    [r,c] = size(theta);
    if ~(r==3 && c==1) && ~(r==1 && c==3)
        if any([r,c]==3)
            if r == 3
                theta = theta(:,1);
            elseif c == 3
                theta = theta(1,:);
            else
                error('weird error')
            end
        else
            error('Function: FullLeg.att_pos_on_demand: input theta vector is not 1x3')
        end
    end

    outPos = zeros(3,1);
    axesMat = zeros(3,3); 
    bodyNum = muscAtt{1,3};

    for i = 1:3
%                 axesMat(:,i) = obj.CR_bodies(:,:,i+1)*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
        axesMat(:,i) = obj.body_obj{i+1}.CR*obj.three_axis_rotation(obj.joint_obj{i}.euler_angs)*[-1;0;0];
    end

    a = axis_angle_rotation(obj,theta(1),axesMat(:,1));
    a = axisAngleRotation(theta(1),axesMat(:,1));
    b = axis_angle_rotation(obj,theta(2),axesMat(:,2));
    c = axis_angle_rotation(obj,theta(3),axesMat(:,3));

    pelPos = obj.organism_position';

    femurpos = pelPos+(obj.body_obj{1}.CR*obj.body_obj{2}.position');
    hiprel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR;

    tibiapos = femurpos+hiprel*obj.body_obj{3}.position';
    kneerel = obj.body_obj{1}.CR*a*obj.body_obj{2}.CR*b*obj.body_obj{3}.CR;
    %kneepos = tibiapos+(kneerel*obj.joint_obj{2}.init_pos);

    footpos = tibiapos+(kneerel*obj.body_obj{4}.position');
    anklerel = kneerel*c*obj.body_obj{4}.CR;
    %anklepos = footpos+(anklerel*obj.joint_obj{3}.init_pos);

    switch bodyNum
        case 1
            outPos = pelPos+obj.body_obj{1}.CR*muscAtt{1};
        case 2
            outPos = (femurpos+hiprel*muscAtt{1});
        case 3
            outPos = (tibiapos+kneerel*muscAtt{1});
        case 4
            outPos = footpos+(anklerel*muscAtt{1});
    end

end