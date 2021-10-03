function [outI,fVal,outFull] = legPositionOpt(desPos,prevPos)
    simPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim";
    trialPath = 'C:\Users\fry16\OneDrive\Desktop\stimulusOptimization\';
    jsystem(['del /q ',trialPath,'*']);

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
    
%     options = optimoptions('fmincon','Algorithm','sqp','TolFun',1e-10,'FiniteDifferenceStepSize',.01,'TolCon',1e-8,'Display','iter-detailed');
     options = optimoptions('fmincon','Algorithm','sqp','Display','iter-detailed',...
         'PlotFcn',{'optimplotx';'optimplotfval';'optimplotstepsize'},'UseParallel',true,'FiniteDifferenceType','central');
    %options = optimoptions('fmincon','Algorithm','sqp','FiniteDifferenceType','central','UseParallel',true,'Display','iter-detailed');
     pattOpts = optimoptions('patternsearch',...
         'PlotFcn',{'psplotbestx','psplotbestf','psplotmeshsize'},'UseParallel',true,'MaxIterations',200,'MeshTolerance',8e-4,'FunctionTolerance',1e-3,'MaxFunEvals',150);
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
        i0 = randi(20,[1 7]);
        %i0 = zeros([1 7]);
    end
    i0 = prevPos;
    i0(i0>=20) = 19.999;
    i0(i0<=0) =.001;
    fun = @(I) objFun(I,inText,stimInds,trialPath,desPos,i0);
    lb = zeros(1,length(stimInds));
    ub = 20.*ones(1,length(stimInds));
    %i0 = [3 14 3 7.75 .7 4 3];
    %i0 = [3.3243   15.1194   18.1518    7.5710    2.3926    5.3407    3.1727];
    %i0 = zeros([1 7]);
    [outI,fVal] = patternsearch(fun,i0,[],[],[],[],lb,ub,[],pattOpts);
    %[outI,fVal] = fmincon(fun,i0,[],[],[],[],lb,ub,[],options);
    outFull = cell(length(stimNames),2);
    for pp = 1:length(stimNames)
        outFull{pp,1} = stimNames{pp};
        outFull{pp,2} = outI(pp);
    end
    [~,outPos] = fun(outI);
    outFull{length(stimNames)+1,1} = outPos(1);
    outFull{length(stimNames)+1,2} = outPos(2);
    outFull{length(stimNames)+1,3} = outPos(3);
    
    function [outFunc,simPos] = objFun(I,inText,stimInds,trialPath,desPos,prevPos)
        %desired position. Make sure the ORDER corresponds to jointProfile order
            %desPos = [-.562388, -.15, .1583246];
            %desPos = [-.487339, -.1888, -.9496332];
        %write V to sim file
            for ii = 1:length(stimInds)
                inText{stimInds(ii)+9} = replaceBetween(inText{stimInds(ii)+9},'>','<',[num2str(I(ii)/10),'e-008']);
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
            executable = [string([sour_folder,'\AnimatSimulator']),string(trialFileName)];
            %executable = ['"',sour_folder,'\AnimatSimulator','" "',trialFileName,'"'];
            %[~,~,err] = jsystem(executable,'noshell');
            jsystem(executable);
        %subtract simPos from desPos
                ds = importdata([trialPath,['JointMotion',['_',num2str(numJM)],'.txt']]);
                jointProfile = ds.data(:,3:5);
                temp = length(jointProfile);
                lbb = floor(temp/3);
                ubb = floor(.99*temp);
                simPos = mean(jointProfile(lbb:ubb,:));
                outFunc = (1/length(simPos)).*sum((simPos-desPos).^2)+(.5/length(I))*sum(((I-prevPos)./(prevPos)).^2);
    end
end