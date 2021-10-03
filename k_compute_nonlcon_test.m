function [c,ceq] = k_compute_nonlcon_test(k,inText,muscDat)
    muscDat(:,3) = num2cell(k);
        % Update the K values
    for ii = 1:length(muscDat)
        muscInd = find(contains(inText,['<Name>',muscDat{ii,1},'</Name>'],'IgnoreCase',true));
        ksInd = find(contains(inText(muscInd:end),'Kse'),1,'first')+muscInd-1;
        kpInd = find(contains(inText(muscInd:end),'Kpe'),1,'first')+muscInd-1;
        factor = muscDat{ii,3}/muscDat{ii,2}; 
        if ~(factor > .99 && factor < 1.01)
            origKs = double(extractBetween(string(inText{ksInd}),'>','</'));
            origKp = double(extractBetween(string(inText{kpInd}),'>','</'));
            inText{ksInd} = replaceBetween(inText{ksInd},'>','</',num2str(origKs*factor));
            inText{kpInd} = replaceBetween(inText{kpInd},'>','</',num2str(origKp*factor));
        end
    end
    
    % Turn off the constant joint stimuli
    cStimNames = {'ConstantHip';'ConstantKnee';'ConstantxAnkle'};
    for ii = 1:length(cStimNames)
        cStimInd = find(contains(inText,cStimNames{ii}));
        inText{cStimInd+2} = '<Enabled>False</Enabled>';
    end
    
    testSim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone_neutral_test.asim";
    fileID = fopen(testSim,'w');
    fprintf(fileID,'%s\n',inText{:});
    fclose(fileID);

    % NW neutral is [103.58, 116.73, 138.11] in NW coords
    % In [5.1427, 14.504, 21.863] in Animatlab coords [0.089757, 0.25314, 0.38158]
    sdata = processSimData(testSim);
    testJM(:,1) = sdata(6).data(:,contains(sdata(6).data_cols,'Hip'));
    testJM(:,2) = sdata(6).data(:,contains(sdata(6).data_cols,'Knee'));
    testJM(:,3) = sdata(6).data(:,contains(sdata(6).data_cols,'Ankle'));
    
    ceq = (mean(testJM(9e3:18e3,:))-[0.089757, 0.25314, 0.38158]).^2;
    c = [];
end