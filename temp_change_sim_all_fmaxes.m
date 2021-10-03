    inprojpath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl.aproj";
    insimpath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_walking.asim";
    intext = importdata(inprojpath);

    fmax_inds = find(contains(intext,'<MaximumTension'));

    for ii = 1:length(fmax_inds)
        fmax_new = 2*double(extractBetween(string(intext{fmax_inds(ii)}),'Actual="','"'));
        intext{fmax_inds(ii)} = replaceBetween(intext{fmax_inds(ii)},'Actual="','"',num2str(fmax_new));
        intext{fmax_inds(ii)} = replaceBetween(intext{fmax_inds(ii)},'Value="','"',num2str(fmax_new));
    end

    % Write simText to an ASIM document
        fileID = fopen(inprojpath,'w');
        fprintf(fileID,'%s\n',intext{:});
        fclose(fileID);
%     
%     % Run the ASIM file    
%         obj = design_synergy(inSimPath);
%     
%     % Find simulation torque
%         [tSraw,~,pmTorque] = compute_passive_joint_torque(obj);