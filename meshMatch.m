function out_text = meshMatch(docpath)
    % Updates all mesh file paths for a given document and overwrites it.
    % Input: docpath: filepath to a Simulation file (.asim) or a Project file (.aproj)
    
    if iscell(docpath)
        doc_text = docpath;
        write2file = 0;
    else
        if ~isfile(docpath)
            error('meshMatch: File not found.')
        end
        %% Store the document text as a cell array
        doc_text = importdata(docpath);
        write2file = 1;
    end
    
    %% Update all MeshFiles
    meshInds = find(contains(doc_text,'<MeshFile>'));
    meshLines = extractBetween(doc_text(meshInds),'<MeshFile>','</MeshFile>');
    
    for ii = 1:length(meshLines)
        line = meshLines{ii};
        slashInds = strfind(line,'\');
        fileTrailer = line(slashInds(end-1):end);
        outline = ['<MeshFile>',pwd,'\Animatlab\Meshes\Rat',fileTrailer,'</MeshFile>'];
        doc_text{meshInds(ii)} = outline;
    end
    
    %% Update all ConvexMeshFiles
    convexInds = find(contains(doc_text,'<ConvexMeshFile>'));
    convexLines = extractBetween(doc_text(convexInds),'<ConvexMeshFile>','</ConvexMeshFile>');
    
    for ii = 1:length(convexLines)
        line = convexLines{ii};
        slashInds = strfind(line,'\');
        fileTrailer = line(slashInds(end-1):end);
        outline = ['<ConvexMeshFile>',pwd,'\Animatlab\Meshes\Rat',fileTrailer,'</ConvexMeshFile>'];
        doc_text{convexInds(ii)} = outline;
    end
    
    %% Write updated document text to the original document
    if write2file == 1
        fid = fopen(docpath,'w');
        fprintf(fid,'%s\n',doc_text{:});
        fclose(fid);
        disp('Meshes updated in file.')
    else
        out_text = doc_text;
    end
end
