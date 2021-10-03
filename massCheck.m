function massCheck(docpath)
    % Sets all mass and density values of bones in a document (.asim or .aproj) to user-defined values.
    % Animatlab will sometimes miscalculates the volume of meshes and that value can't be changed.
    % To make the volume make sense, it automatically sets the scale of the density and mass.
    % Input: docpath: filepath to a Simulation file (.asim) or a Project file (.aproj)
    
    if ~isfile(docpath)
        error('meshMatch: File not found.')
    else
        if contains(docpath,'.aproj')
            type = 'proj';
        elseif contains(docpath,'.asim')
            type = 'sim';
        else
            error('Input file is not .asim or .aproj')
        end
    end
    
    %% Store the document text as a cell array
    doc_text = importdata(docpath);
    % Find the indices of the bones
    boneInds = find(contains(doc_text,{'Pelvis</Name>','Femur</Name>','Tibia</Name>','Foot</Name>'}));
    
    % The correct values of the bone masses in (g)
    correctMass = [31.824;...
                   14.141;...
                   3.339;...
                   1.571];
   
   % The correct values of the bone densities in (g/cm^3)
    correctDens = [3;...
                   8.26;...
                   2.471;...
                   3];
    
    if isempty(boneInds)
        error('Cannot locate bone indices in document')
    end
    
    switch type
        case 'proj'
            massBound = '<Mass Value=';
            densityBound = '<Density Value=';
        case 'sim'
            massBound = '<Mass>';
            densityBound = '<Density>';
    end
    
    for ii = 1:size(boneInds,1)
        bI = boneInds(ii);
        mInd = find(contains(doc_text(bI:end),massBound),1,'first')+bI-1;
        dInd = find(contains(doc_text(bI:end),densityBound),1,'first')+bI-1;
        switch type
            case 'proj'
                doc_text{mInd} = ['<Mass Value="',num2str(correctMass(ii)),'" Scale="None" Actual="',num2str(correctMass(ii)),'"/>'];
                doc_text{dInd} = ['<Density Value="',num2str(correctDens(ii)),'" Scale="None" Actual="',num2str(correctDens(ii)),'"/>'];
            case 'sim'
                doc_text{mInd} = ['<Mass>',num2str(correctMass(ii)),'</Mass>'];
                doc_text{dInd} = ['<Density>',num2str(correctDens(ii)),'</Density>'];
        end
    end
    
    %% Write updated document text to the original document
    fid = fopen(docpath,'w');
    formatSpec = '%s\n';
    nrows = size(doc_text,1);
    for row = 1:nrows
        fprintf(fid,formatSpec,doc_text{row,:});
    end
    fclose(fid);
end
