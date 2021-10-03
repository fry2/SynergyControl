function [objCell,objTorques_act,objTorques_tot] = compare_rat_and_cat_walking()
    fileString = @(size,walk) ['G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_Size_',size,'_Walk_',walk,'.asim'];
    oneHundredizer = @(inMat,trans) [interp1(1:length(inMat(1:trans,:)),inMat(1:trans,:),linspace(1,length(inMat(1:trans,:)),37)),...
                                     interp1(1:length(inMat(trans:end,:)),inMat(trans:end,:),linspace(1,length(inMat(trans:end,:)),63))];
    oneHundredizer_mat = @(inMat,trans) [interp1(1:length(inMat(1:trans,:)),inMat(1:trans,:),linspace(1,length(inMat(1:trans,:)),37));...
                                     interp1(1:length(inMat(trans:end,:)),inMat(trans:end,:),linspace(1,length(inMat(trans:end,:)),63))];
    dataNames = {'Rat';'Cat';'Horse'}; counter = 1; objCell = cell(1,length(dataNames)^2); objTorques = objCell;
    for ii = 1:length(dataNames)
        for jj = 1:length(dataNames)
            fPath = fileString(dataNames{ii},dataNames{jj});
            objCell{counter} = design_synergy(fPath);
            results_cell = pedotti_optimization(objCell{counter});
            force_mn = results_cell{3,2};
            if strcmp(dataNames{jj},'Cat')
                beg = 385; mid = 724; ennd = 1056; % is cat
            elseif strcmp(dataNames{jj},'Horse')
                beg = 895; mid = 1479; ennd = 2300; % is horse
            else
                beg = 654; mid = 1034; ennd = 1676;  % is else (rat)
            end            
            rm = find_relevant_muscles(objCell{counter});
            for kk = 1:3
                [moment_arms] = compute_joint_moment_arms(objCell{counter},kk,1);
                moment_arms = moment_arms(rm{kk},:)'./1000;
                active_forces = force_mn(:,rm{kk});
                torque_temp = oneHundredizer(sum(moment_arms(beg:ennd,:).*active_forces(beg:ennd,:),2),mid-beg);
%                 t_swing = -1 + 2.*(torque_temp(1:37) - min(torque_temp(1:37)))./(max(torque_temp(1:37)) - min(torque_temp(1:37)));
%                 t_stance = -1 + 2.*(torque_temp(38:end) - min(torque_temp(38:end)))./(max(torque_temp(38:end)) - min(torque_temp(38:end)));
                t_swing = torque_temp(1:37)./max(abs(torque_temp(1:37))); t_stance = torque_temp(37:end)./max(abs(torque_temp(37:end)));
                objTorques_act{counter}(:,kk) = [t_swing,t_stance];
            end
            tjt = compute_total_joint_torque(objCell{counter},0);
            tjt = oneHundredizer_mat(tjt(beg:ennd,:),mid-beg);
%             tjt_swing = -1 + 2.*(tjt(1:37,:) - min(tjt(1:37,:)))./(max(tjt(1:37,:)) - min(tjt(1:37,:)));
%             tjt_stance = -1 + 2.*(tjt(38:end,:) - min(tjt(38:end,:)))./(max(tjt(38:end,:)) - min(tjt(38:end,:)));
            tjt_swing = tjt(1:37,:)./max(abs(tjt(1:37,:))); tjt_stance = tjt(37:end,:)./max(abs(tjt(37:end,:)));
            objTorques_tot{counter} = -1.*[tjt_swing;tjt_stance];
            disp(['Finished obj ',num2str(counter),' of ',num2str(length(dataNames)^2),' - ',dataNames{ii},' walking like a ',dataNames{jj}]);
            counter = counter + 1;
        end
    end
end