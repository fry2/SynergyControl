function [outVal,outObj] = temp_find_servogain(inVec,inText,goodWave)
    mvInds = find(contains(inText,'<MaxVelocity>'));
    fmInds = mvInds-1;
    sgInds = find(contains(inText,'<ServoGain>'));
    
    for ii = 1:length(fmInds)
        inText{fmInds(ii)} = replaceBetween(inText{fmInds(ii)},'>','</',num2str(inVec(1)));
    end
    
    for ii = 1:length(mvInds)
        inText{mvInds(ii)} = replaceBetween(inText{mvInds(ii)},'>','</',num2str(inVec(2)));
    end
    
    for ii = 1:length(sgInds)
        inText{sgInds(ii)} = replaceBetween(inText{sgInds(ii)},'>','</',num2str(inVec(3)));
    end
    
    % Write simText to an ASIM document
    jobSavePath = 'C:\Users\fry16\OneDrive\Documents\JointDampingOpt\InjectedProject\JointDampingOpt_injected_Standalone_servoGain.asim';
    fileID = fopen(jobSavePath,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);
    
    outObj = design_synergy(jobSavePath);
    outVal = sum(sum((outObj.theta_motion-goodWave).^2));
end