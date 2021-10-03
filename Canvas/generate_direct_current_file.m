function [outData,stimPath] = generate_direct_current_file(inObj,inWave,stimPath)
    % For an input data vector, generate a .txt file that Animatlab can read from using a
    % inverse muscle dynamics current stimulus

    % Input (inObj) struct: a CanvasModel object (ex. nsys)
    % Input (inWave) Nx1 double array: vector of current stimulus waveform
    % Input (stimName) char: character name of stimulus to the waveform to
    % Output (outData) Nx2 double array: two column array with a time and data column
    % Output (stimPath) char: character array of the .txt file data path

    inWave(isnan(inWave)) = 0;
    [r,c] = size(inWave);
    if c>r
        inWave = inWave';
    end

    % The input current waveform must be in units of A. Typically we deal in units of [0,20] nA, so it is important to make this change.
    if max(inWave)>20e-9    
        inWave = inWave*1e-9;
    end

    dt = inObj.proj_params.physicstimestep*1e-3;
    tEnd = inObj.proj_params.simendtime;
    timeVec = 0:dt:tEnd;
    
    % we need to upsample the input wave to match the time vector
    inWave = interp1(1:length(inWave),inWave,linspace(1,length(inWave),length(timeVec)));
    
    outData = [timeVec',inWave'];
    
    outText = cell(length(inWave)+1,1);
    outText{1} = ['Time',sprintf('\t'),'Current'];
    
    strTime = cellstr(num2str(timeVec','%.5f\t\n'));
    strCurrent = cellstr(num2str(inWave','%.4e\t\n'));

    tCell = cell(length(timeVec),1);
    [tCell{:}] = deal(sprintf('\t'));
    outText(2:length(timeVec)+1,1) = strcat(strcat(strTime,tCell),strCurrent);
    
%     strTime = cellfun(@strtrim,cellstr(string(num2str(timeVec'))),'UniformOutput',false);
%     strCurrent = cellfun(@strtrim,cellstr(string(num2str(inWave'))),'UniformOutput',false);
%     outText(2:length(strTime)+1,:) = [strTime,strCurrent];
%     
%     for ii1 = 2:length(strTime)+1
%         outText{ii1,1} = [strTime{ii1-1},sprintf('\t'),strCurrent{ii1-1}];
%     end

    stimPath = [pwd,'\Data\InjectedCurrent\',stimPath,'.txt'];
    fileID = fopen(stimPath,'w');
    fprintf(fileID,'%s\n',outText{:});
    fclose(fileID);
end