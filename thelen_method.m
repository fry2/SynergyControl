% Following Thelen 2003 "Generating dynamic simulations of..."
curvars = whos;
if ~any(contains({curvars.name},'obj'))
    [obj,sim_file,joints,bodies,joint_limits,joint_profile,sdata] = design_synergy("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim");
end
%% Stage 1: Desired Accelerations
startInd = 120;
% time = joint_profile(startInd:end-10,1);
% q0_exp = joint_profile(startInd:end-10,2:end);
time = obj.theta_motion_time(startInd:end);
q0_exp = obj.theta_motion(startInd:end,:);

[~,q1_exp] = gradient(q0_exp,obj.dt_motion);
samplingFreq = 1/obj.dt_motion;
[b,a] = butter(4,6/(samplingFreq/2),'low');
q1_exp = filtfilt(b,a,q1_exp);

[~,q2_exp] = gradient(q1_exp,obj.dt_motion);
samplingFreq = 1/obj.dt_motion;
[b,a] = butter(4,6/(samplingFreq/2),'low');
q2_exp = filtfilt(b,a,q2_exp);

kv = 40;
kp = 400;
ku = 10;

%% Stage 2
    %[Mout,G,C] = MLS_maker(obj);
    % pattOpts = optimoptions('patternsearch',...
    %         'PlotFcn',{'psplotbestx','psplotbestf','psplotmeshsize'},'UseParallel',true,'FunctionTolerance',1e-4,'MeshTolerance',1e-4);
               pattOpts = optimoptions('patternsearch','Display','iter','UseParallel',true,'FunctionTolerance',1e-4,'MeshTolerance',1e-3);
    % Prepare a matrix fo muscle information for converting steady state force to steady state activation
    numMuscles = length(obj.musc_obj);
    musc_info = zeros(numMuscles,3);
    for ii = 1:numMuscles
        muscle = obj.musc_obj{ii};
        musc_info(ii,1) = muscle.Kse;
        musc_info(ii,2) = muscle.Kpe;
        musc_info(ii,3) = muscle.RestingLength;
        musc_info(ii,4) = muscle.max_force;
        musc_info(ii,5) = muscle.ST_max;
        musc_info(ii,6) = muscle.steepness;
        musc_info(ii,7) = muscle.x_off;
    end
    
    q2_exp = q2_exp(4374:10560,:);
    samplingInds = 1:length(q2_exp);
    samplingInds = floor(linspace(1,length(q2_exp),30));
    %samplingInds = 1:20;
    
    %%
    clear I fs
    counter = 1; 
for ii = samplingInds
    % Stage 1
    if ii == 1
        q0_feedback = q0_exp(ii,:)';
        q1_feedback = q1_exp(ii,:)';
    end
    q2_d = q2_exp(ii,:)'+kv*(q1_exp(ii,:)'-q1_feedback)+kp*(q0_exp(ii,:)'-q0_feedback);
    
    % Stage 2: Optimization Process
    %fun = @(fs) thelen_method_obj_func(fs,musc_info);
    fun = @(fs) sum((fs./musc_info(:,4)).^2);
    lb = zeros(numMuscles,1);
    for jj = 1:38
        Ks = musc_info(jj,1); Kp = musc_info(jj,2); Lr = musc_info(jj,3); STmax = musc_info(jj,5);
        ub(jj) = (Ks/(Ks+Kp)).*(1+(Kp.*Lr./(4.*STmax)).^2).*STmax;
        ub(jj) = STmax;
    end
    %ub = musc_info(:,4);
%     if ii == 1
%         fs0 = zeros(numMuscles,1);
%     else
%         fs0 = fs(:,counter-1);
%     end
    
%     [Mout,Mprime] = compute_mass_matrix(obj,q0_exp(ii,:));
%     Mvec(:,counter) = Mout*q2_exp(ii,:)';
%     G = compute_gravity_vector(obj,q0_exp(ii,:))';
%     Gvec(:,counter) = G;
%     C = compute_coriolis_matrix(obj,q0_exp(ii,:),q1_exp(ii,:));
%     Cvec(:,counter) = C*q1_exp(ii,:)'.^2;
%     R = leg_moment_arms(obj,q0_exp(ii,:))'./1000;
%     test1(counter) = R(1,1);
%     
%     beq(:,counter) = -(Mout*q2_exp(ii,:)'-G-C*q1_exp(ii,:)'.^2);
%     Aeq = R;
    %[fs(:,counter),fVal] = patternsearch(fun,fs0,[],[],Aeq,beq(:,counter),lb,ub,[],pattOpts);
    %[fs(:,counter)] = fmincon(fun,fs0,[],[],Aeq,beq(:,counter),lb,ub,[],optimoptions('fmincon','Display','none','Algorithm','sqp','OptimalityTolerance',1e-4));
    
    outLens = muscle_lengths_on_demand(obj,q0_exp(ii,:));
    M = compute_mass_matrix(obj,q0_exp(ii,:));
    G = compute_gravity_vector(obj,q0_exp(ii,:));
    R = leg_moment_arms(obj,q0_exp(ii,:))'./1000;
    fun2 = @(Am) qdot_dot(obj,Am,q2_exp(ii,:),outLens,M,G,R);
    if ii == 1
        Am0 = zeros(numMuscles,1);
    else
        Am0 = Am_out(:,counter-1);
    end
    %[Am_out(:,counter)] = fmincon(fun2,Am0,[],[],[],[],lb,ub,[],optimoptions('fmincon','Display','none','Algorithm','sqp','OptimalityTolerance',1e-4));
     temp = patternsearch(fun2,Am0,[],[],[],[],lb,ub,[],pattOpts);
     Am_out(:,counter) = temp;
%     fs(fs(:,counter)<1e-5,counter) = 0;
%     [outVal,A] = thelen_method_obj_func(fs(:,counter),musc_info);
%     
%     % Stage 3: Find neuron stimuli to generate activations
%     if ii == 1
%         A_feedback = A;
%     end
%     %Ainput = A + ku*(A-A_feedback);
%     I(:,counter) = thelen_method_current_solve(A,musc_info);
    disp(num2str(counter))
    counter  = counter +1;  
    %keyboard
    
    % Stage 4
    %q2_feedback = inv(Mout(:,:,ii))*(G(ii,:)'+C(:,:,ii)+Aeq*fs);
end

function outQ = qdot_dot(obj,Am,q2des,L,M,G,R)
    Tens = @(L,Am,Ks,Kp,Lr) (Ks/(Ks+Kp))*(max(Kp*(L-Lr),0)+(1-4*(L-Lr).^2/(Lr^2))*Am);
    numMusc = length(obj.musc_obj);
    fm = zeros(numMusc,1);
    %Am = fm;
    for ii = 1:numMusc
        muscle = obj.musc_obj{ii};
        Lr = muscle.RestingLength; Kp = muscle.Kpe; Ks = muscle.Kse;
        fm(ii) = Tens(L(ii),Am(ii),Ks,Kp,Lr);
    end
    outQ = sum((inv(M)*(G'+R*fm)-q2des').^2);
end