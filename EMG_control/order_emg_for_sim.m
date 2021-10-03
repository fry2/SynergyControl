function outEMGdata = order_emg_for_sim(inSimPath,inEMGdata)
    % Function: emg data and simulation files do not always have muscles ordered the same way. To make life easier, this function preprocesses the data
    % so that the data and file order matches. This way, we don't need to do any mapping or worry about pointing to the wrong data
    
    % Input: inSimPath (char or string): path to a simulation file
    % Input: inEMGdata (cell): 2x1 or 1x2 cell that contains EMG waveforms (double) in cell 1 and column names in cell 2
    % Output: outEMGdata (cell): cell of same size as inEMGdata that contains ordered data

    if iscell(inEMGdata)
        [r,c] = size(inEMGdata);
        if (r==2 && c==1) || (r==1 && c==2)
            emgMuscNames = inEMGdata{2};
        else
            error('Input inEMGdata must have dimensions 1x2 or 2x1.')
        end
    else
        error('Input inEMGdata must be a cell.')
    end
        
    %does the inSimPath file exist?
        %does the inSimPath file contain muscle objects?
    if isfile(inSimPath)
        muscle_info = scrapeFileForMuscleInfo(inSimPath);
    else
        error('Input inSimPath does not point to an existing .asim file.')
    end
    
    if length(muscle_info) ~= length(emgMuscNames)
        warning('Warning: order_emg_for_sim: input emg data and input simulation file appear to have different muscles.')
    end
    
    nameKey =   {'VL'  ,    'lh_vastuslateralis'      ;...
                'TA'  ,    'lh_tibialisanterior'      ;...
                'GA'  ,    'lh_gracilisanterior'      ;...
                'RF'  ,    'lh_rectusfemoris'         ;...
                'BFpR',    'lh_bicepsfemorisposterior';...
                'VM'  ,    'lh_vastusmedialis'        ;...
                'SM'  ,    'lh_semimembranosus'       ;...
                'IL'  ,    'lh_illiopsoas'            ;...
                'GS'  ,    'lh_gemellussuperior'      ;...
                'VI'  ,    'lh_vastusintermedius'     ;...
                'GRr' ,    'lh_gracilisposterior'     ;...
                'ST'  ,    'lh_semitendinosus'        };
            
    outEMGdata = cell(size(inEMGdata));
    for ii = 1:length(muscle_info)
        muscName = muscle_info{ii,2};
        if any(contains(nameKey,muscName),'all')
            [mar,~] = find(contains(nameKey,muscName));
            muscAcro = nameKey{mar,1};
            muscName = nameKey{mar,2};
        end
        if any(contains(emgMuscNames,muscAcro),'all') || any(contains(emgMuscNames,muscName),'all')
            [~,emgc] = find(contains(emgMuscNames,muscAcro));
            if isempty(emgc)
                [~,emgc] = find(contains(emgMuscNames,muscName));
            end
        end
        outEMGdata{1}(:,ii) = inEMGdata{1}(:,emgc);
        outEMGdata{2}{ii} = muscName;
    end
end