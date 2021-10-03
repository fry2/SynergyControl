function sdata = temp_stims2sdata(simContents,tflStim,abStim)
        simPath = [pwd,'\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim'];
        tflInd = find(contains(simContents,'stTC1-1TFL'))+13;
        abInd = find(contains(simContents,'stTC1-7AB'))+13;
        simContents{tflInd} = replaceBetween(simContents{tflInd},'>','<',[num2str(tflStim/10),'e-008']);
        simContents{abInd} = replaceBetween(simContents{abInd},'>','<',[num2str(abStim/10),'e-008']);
        fileID = fopen(simPath,'w');
        fprintf(fileID,'%s\n',simContents{:});
        fclose(fileID);
        sdata = processSimData(simPath);
end