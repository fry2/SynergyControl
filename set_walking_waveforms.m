function set_walking_waveforms(inPath)
    % This function will update a given ASIM or APROJ file with a specific, saved walking profile. Additional profiles can be added.
    if contains(inPath,'.asim')
        is_sim = 1;
    elseif contains(inPath,'.aproj')
        is_sim = 0;
    else
        error('Input file path does not point to an ASIM or APROJ file')
    end
    inText = importdata(inPath);
    outText = inText;
    hipInd = find(contains(inText,'<Name>Walking_Hip</Name>'));
    kneeInd = find(contains(inText,'<Name>Walking_Knee</Name>'));
    ankleInd = find(contains(inText,'<Name>Walking_xAnkle</Name>'));
%     
%     switch 1
%         case 1
%             % This option is for normal walking, as defined by Fischer's 02 paper
%             load([pwd,'\Data\walkingWaveforms_normal.mat'],'normalequations','normalequations_proj')
%             simeqs = normalequations; projeqs = normalequations_proj;
%         case 2
%             % This option uses the same data as normal walking, but initializes all joint profiles to start at 0
%             load([pwd,'\Data\walkingWaveforms_init0.mat'],'init0equations','init0equations_proj')
%             simeqs = init0equations; projeqs = init0equations_proj;
%         otherwise
%             error('Choose valid case')
%     end
    
    load([pwd,'\Data\walking_waveforms.mat'],'walking_waveforms')
    simeqs = walking_waveforms.Cat_walking(:,2);
    projeqs = walking_waveforms.Cat_walking(:,1);

    if is_sim
        indmod = 10; eqs = simeqs;
    else
        indmod = 8; eqs = projeqs;
    end
    
    outText{hipInd+indmod} = replaceBetween(inText{hipInd+indmod},'>','</',eqs{1});
    outText{kneeInd+indmod} = replaceBetween(inText{kneeInd+indmod},'>','</',eqs{2});
    outText{ankleInd+indmod} = replaceBetween(inText{ankleInd+indmod},'>','</',eqs{3});
    
    %savePath = inPath;
    savePath = [char(a),'\',char(b),'_cat',char(c)];
    fileID = fopen(savePath,'w');
    fprintf(fileID,'%s\n',outText{:});
    fclose(fileID);
end