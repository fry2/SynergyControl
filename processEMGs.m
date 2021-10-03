function outEMGs = processEMGs
    % Interpolate given EMG signals for a desired sampling rate
    % Input: inEMGs: struct with individual muscle waveforms over a [0,100]% gait profile
    
        % Data from Dr. Matt Tresch Northwestern student thesis Hsin-Yun Yeh "EMG Analysis of Rat Hindlimb in Different Context:
            ... from Downslope to Upslope Walking" June 2019, p.20
    inEMGs = load([pwd,'\Data\Tresch_EMGs.mat']);
    
    EMGname = fieldnames(inEMGs);
    outEMGs = struct();
    desiredSamples = 100;
    for ii = 1:length(EMGname)
        inSig = inEMGs.(EMGname{ii});
        n = length(inSig);
        outSig = interp1(1:n,inSig(:,2),linspace(1,n,desiredSamples))';
        outSig = smoothdata(outSig,'gaussian',.1*desiredSamples);
        outEMGs.(EMGname{ii}) = outSig;
    end
end