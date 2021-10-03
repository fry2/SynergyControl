function [outk,r2scores,recompiled,W,H] = NMFsyncounter(forces)
    % Determine the minimum number of synergies to meet a matching criteria between the recombined and input waveforms.
    % Input: forces: array of forces over time
    % Output: outk: integer number of synergies necessary
    % Output: r2scores: (nx1) array of r^2 values representing similarity of recombined waveforms to input waveforms
    % Output: recompiled: array of size(forces) representing forces as generated from synergy components
    % Output: W: spatial array of synergies representing muscle weighting
    % Output: H: timing array of synergies representing activation over time
    
    if size(forces,1)>size(forces,2)
        forces = forces';
    end
    
    r2scores = zeros(38,1);
    outk = 0;
    
    while any(r2scores < .5) && outk < 25
        outk = outk+1;
        [r2scores,recompiled,W,H] = NMFdecomposition(outk,forces,0);
    end
    
    if outk <= 10
        plotWH(forces,W,H,0);
    else
        disp('Over 10 synergies')
    end
end