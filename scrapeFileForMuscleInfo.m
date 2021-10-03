function muscle_info = scrapeFileForMuscleInfo(fPath)
    % Scrape text for muscle information
    % Input (fPath) char: path to .asim or .aproj file
    % Output: muscle_out Nx3: cell array containing columns of muscle indices, muscle names, and muscle ID's
    
    % If given a simPath to a file that doesn't exist, check if a project exists by the same name and extract muscle info from that
    if ~isfile(fPath)
        if ~ischar(fPath)
            fPath = char(fPath);
        end
        [docPath,docName] = fileparts(fPath);
        if contains(docName,'Standalone')
            temp = strfind(docName,'_Standalone');
            docName = docName(1:temp-1);
        end
        fPath = [docPath,'\',docName,'.aproj'];
    end

    original_text = importdata(fPath);
    muscleInds = find(contains(original_text,'<PartType>AnimatGUI.DataObjects.Physical.Bodies.LinearHillMuscle</PartType>'))-2;
    if isempty(muscleInds)
        error('No muscles in input text')
    end
    muscle_info = cell(length(muscleInds),3);
    for ii = 1:length(muscleInds)
        muscle_info{ii,1} = muscleInds(ii);
            mNametrim = strtrim(original_text{muscleInds(ii)-1});
        muscle_info{ii,2} = lower(mNametrim(7:end-7));
            mIDtrim = strtrim(original_text{muscleInds(ii)});
        muscle_info{ii,3} = mIDtrim(5:end-5);
    end
end