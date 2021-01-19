function [newTens,origTens,newJM,origJM] = paramChanger(inSimPath,scales)
    % For a given simulation file, scale all muscle muscle parameters by a certain value
    % inSimPath: char: simulation path
    % scales: double (1x3): the amount to scale [B, Ks, Kp] % ex: [1 1 10]
    sdata = processSimData(inSimPath);
    tensionInd = contains({sdata.name},'KeyMuscTen');
    origTens = sdata(tensionInd).data;
    
    jmInd = contains({sdata.name},'JointMotion');
    origJM = sdata(jmInd).data.*(180/pi)+[98.4373 102.226 116.2473];
    
    origtext = importdata(inSimPath);
    muscInds = find(contains(origtext,'<Type>LinearHillMuscle</Type>'));
    for ii = 1:length(muscInds)
        % B
        oldPar = double(extractBetween(string(origtext{muscInds(ii)+45}),'>','<'));
        origtext{muscInds(ii)+45} = replaceBetween(origtext{muscInds(ii)+45},'>','<',num2str(oldPar*scales(1)));
        % Ks
        oldPar = double(extractBetween(string(origtext{muscInds(ii)+43}),'>','<'));
        origtext{muscInds(ii)+43} = replaceBetween(origtext{muscInds(ii)+43},'>','<',num2str(oldPar*scales(2)));
        % Kp
        oldPar = double(extractBetween(string(origtext{muscInds(ii)+44}),'>','<'));
        origtext{muscInds(ii)+44} = replaceBetween(origtext{muscInds(ii)+44},'>','<',num2str(oldPar*scales(3)));
    end
    
    txtInd1 = find(contains(origtext,'OutputFilename>Joint'));
    origtext{txtInd1} = replaceBetween(origtext{txtInd1},'ion','.txt','_temp');
    txtInd2 = find(contains(origtext,'OutputFilename>KeyMuscTen'));
    origtext{txtInd2} = replaceBetween(origtext{txtInd2},'Ten','.txt','_temp');
    
    jobSavePath = [pwd,'\temp1.asim'];
    fileID = fopen(jobSavePath,'w');
    fprintf(fileID,'%s\n',origtext{:});
    fclose(fileID);
    
    newData = processSimData(jobSavePath);
    tensionInd2 = contains({newData.name},'KeyMuscTen');
    newTens = newData(tensionInd2).data;
    
    jmInd = contains({newData.name},'JointMotion');
    newJM = newData(jmInd).data.*(180/pi)+[98.4373 102.226 116.2473];
    
    figure('Position',[962,2,958,994]);
    subplot(4,1,1)
    plot(origJM); title('paramChanger output: Original Joint Motion')
    subplot(4,1,2)
    plot(newJM); title('New Joint Motion (w/ scaled vals)')
    subplot(4,1,3)
    plot(origTens); title('Original Tension')
    subplot(4,1,4)
    plot(newTens); title('New Tension (w/ scaled vals)')
    
    delete(jobSavePath)
end