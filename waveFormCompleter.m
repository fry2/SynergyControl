function completeWaves = waveFormCompleter(inWaves,linearScale,meanScale)
    % Extrapolate missing joint motion data from existing data, like the data in spreadsheets Alex processed in the ParameterCalculation
        % (KinematicDataProcessing)
    % Input: inWaves: (timesteps x trials) raw input waveforms for a single joint, both complete and incomplete [e.g. BackRaw(:,:,1)]
    % Input: linearScale: how adherent the waveform extrapolator is to maintaining the linear projection of the existing data of incomplete waves. In dimensionless 
        %units of "steps" that are contingent on the length of the input data (e.g. 500 for input waveform length of 1000)
    % Input: meanScale: how adherent the waveform is to the mean of existing data. The idea being that extrapolated data will become increasingly similar to
        % mean data as it distances itself from its existing wave data
    % Output: completeWaves: waveforms with dimensions of inWaves except all waveforms are now complete
    inSz = size(inWaves);
    
    if nargin ~= 3
        linearScale = 400;
        meanScale = 400;
    end
    
    if any(isnan(inWaves(1,:)))
        inWaves = flipud(inWaves);
        nanBool = isnan(inWaves(end,:,1));
        flipper = 1;
    else
        nanBool = isnan(inWaves(1,:,1));
    end
    
    completeWaves = zeros(inSz);
    completeWaves(:,~nanBool) = inWaves(:,~nanBool);

     for i = 1:inSz(2)
        if nanBool(i)
            startInd = find(isnan(inWaves(:,i)),1,'first');
            extWave = inWaves(1:startInd-1,i);
            scaleVec = zeros(1,startInd-1);
            linVec = ones(1,startInd-1);
            for j = startInd:inSz(1)
                adder = mean(inWaves(j-1,~isnan(inWaves(j-1,:))))-extWave(j-1);
                adScaler = min(exp((7/meanScale)*(j-startInd)-7),1);
                linearProj = (extWave(startInd-1)-extWave(startInd-2));
                linScaler = max((-1/linearScale)*(j-startInd)+1,0);
                scaleVec(j) = adScaler;
                linVec(j) = linScaler;
                 extWave(j) = extWave(j-1)+adScaler*adder+linScaler*linearProj;
            end
            completeWaves(:,i) = extWave;
        end
    end
    if flipper
        completeWaves = flipud(completeWaves);
        inWaves = flipud(inWaves);
    end
    
    to_plot = 1;
    if to_plot
        subplot(2,1,1)
        plot(inWaves)
        xlabel('inWaves')
        xlim([0 inSz(1)])
        subplot(2,1,2)
        plot(completeWaves)
        xlabel('completeWaves')
        xlim([0 inSz(1)])
    end
end