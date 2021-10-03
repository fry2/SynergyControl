function [outI,fVal,outFull] = legWaveformOpt(desWave)
    simPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim";
    trialPath = 'C:\Users\fry16\OneDrive\Desktop\stimulusOptimization\';
    delete([trialPath,'*.*'])

    inText = importdata(simPath);
    stimInds = find(contains(inText,'<CurrentType>Tonic</CurrentType>'));
    stimNames = inText(stimInds-4);
    %check that all the simuli link to muscles
        for jj = 1:length(stimInds)
            targetStr = extractBetween(string(inText{stimInds(jj)+2}),'>','<');
            targetInd = find(contains(inText,targetStr),1,'first');
            if ~strcmp(inText{targetInd-1},'<Neuron>') || ~strcmp(inText{targetInd+2},'<Enabled>True</Enabled>')
                error('Stimulus Error')
            end
        end
    to_plot = 0;
    fun = @(I) objFun(I,inText,stimInds,trialPath,desWave,to_plot);
    tester = randi(20,[1,6*7]);
    options = optimoptions('fmincon','Algorithm','sqp','UseParallel',true,'Display','iter-detailed',...
         'PlotFcn',{'optimplotx';'optimplotfval';'optimplotstepsize'},'FiniteDifferenceType','forward','FiniteDifferenceStepSize',.5);
    findIO = 0;
    if findIO
        numTests =  20;
        counter = 1;
        bestVal = 1000;
        while counter <= numTests
            %testI = 20.*rand([1 7]);
            testI = [1.1570   18.0857    7.7714    7.5301    0.4140    4.8178    6.2406]+(7.*rand([1 7]) - 3.5);
            testI(testI<0) = 0;
            testOut = fun(testI);
            if testOut < bestVal
                bestVal = testOut;
                i0 = testI;
            end
            counter = counter+1;
            if ~mod(counter,20)
                disp(num2str(counter))
            end
        end
        delete([trialPath,'*.*'])
    else
        i0 = randi(20,[1 6*7]);
        %i0 = zeros([1 7]);
    end
    i0(i0==20) = 19.999;
    i0(i0==0) =.001;
    %lb = zeros(1,length(stimInds));
%     lb = -1000.*ones(1,length(i0));
    lb = repmat([0 0 -Inf],[1 14]);
    ub = 1000.*ones(1,length(i0));
    %i0 = [3 14 3 7.75 .7 4 3];
    %i0 = [3.3243   15.1194   18.1518    7.5710    2.3926    5.3407    3.1727];
    %i0 = zeros([1 7]);
    [outI,fVal] = fmincon(fun,i0,[],[],[],[],lb,ub,[],options);
    outFull = 1;
    objFun(outI,inText,stimInds,trialPath,desWave,1);
    return
    sampInds = round(linspace(1,length(desWave),50));
    outI = zeros(7,length(sampInds));
    i0 = 3*ones([1 7]);
    for ii1 = 1:length(sampInds)
        if ii1 > 1
            i0 = outI(:,ii1-1);
        end
        fun = @(I) objFun(I,inText,stimInds,trialPath,desWave(sampInds(ii1),:));
        [outI(:,ii1),fVal] = fmincon(fun,i0,[],[],[],[],lb,ub,[],options);
        disp([num2str(ii1),num2str(fVal)])
        disp(['[',num2str(ii1),',',num2str(fVal),']'])
    end
        outFull = cell(length(stimNames),2);
        return
    for pp = 1:length(stimNames)
        outFull{pp,1} = stimNames{pp};
        outFull{pp,2} = outI(pp);
    end
    [~,outPos] = fun(outI);
    outFull{length(stimNames)+1,1} = outPos(1);
    outFull{length(stimNames)+1,2} = outPos(2);
    outFull{length(stimNames)+1,3} = outPos(3);
    
    function [outFunc,simPos] = objFun(I,inText,stimInds,trialPath,desPos,to_plot)
        %desired position. Make sure the ORDER corresponds to jointProfile order
            %desPos = [-.562388, -.15, .1583246];
            %desPos = [-.487339, -.1888, -.9496332];
        %write V to sim file
        I = reshape(I,[6 7]);
            coeffs = cell(6,2);
            coeffs(:,1) = {'a1';'b1';'c1';'a2';'b2';'c2'};
            for ii = 1:length(stimInds)
                coeffs(:,2) = num2cell(I(:,ii)');
                sim_eqn = sum_of_sines_maker(coeffs,0);
%                 inText{stimInds(ii)+9} = replaceBetween(inText{stimInds(ii)+9},'>','<',[num2str(I(ii)/10),'e-008']); % To replace a constant
                inText{stimInds(ii)+12} = replaceBetween(inText{stimInds(ii)+12},'>','<',sim_eqn);
            end
        %specify the output JointMotion file name
            txtInd = find(contains(inText,'OutputFilename>Joint'));
            numJM = randi(100000000);
            inText{txtInd} = replaceBetween(inText{txtInd},'ion','.txt',['_',num2str(numJM)]);
        %write content to file
            trialFileName = [trialPath,'optimTest_',num2str(numJM),'.asim'];
            fileID = fopen(trialFileName,'w');
            fprintf(fileID,'%s\n',inText{:});
            fclose(fileID);
        %evaluate the sim file
            sour_folder = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin';
            executable = ['"',sour_folder,'\AnimatSimulator" "',trialFileName,'"'];
            [~,~,err] = jsystem(executable,'noshell');
        %subtract simPos from desPos
            if isempty(err)
                ds = importdata([trialPath,['JointMotion',['_',num2str(numJM)],'.txt']]);
                jointProfile = ds.data(:,3:5);
                %temp = length(jointProfile);
                %lbb = floor(temp/3);
                %ubb = floor(.99*temp);
                %simPos = mean(jointProfile(lbb:ubb,:));
                 minLen = min([length(jointProfile),length(desPos)]);
                 outFunc = sum(sum((desPos(1:minLen,:)-jointProfile(1:minLen,:)).^2),2);
                %outFunc = (simPos(1)-desPos(1))^2+(simPos(2)-desPos(2))^2+(simPos(3)-desPos(3))^2;
            end
        if to_plot
            figure;
            subplot(2,1,1)
                plot(desPos)
            subplot(2,1,2)
                plot(jointProfile)
        end
    end
end