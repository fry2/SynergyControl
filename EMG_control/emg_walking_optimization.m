function [outScale,fVal] = emg_walking_optimization(inEMGdata,jointMotion,initialPoint)
%%
    % Ensure that input EMG data is ordered the same as the simulation file
    simPath = 'G:\My Drive\Rat\SynergyControl\Data\EMG_walking_baseSim.asim';
    % For this to work, you need a .asim file that contains existing DC stimulus' for each muscle
    % The optimization is sped up by having these already established and simply changing the
    % Stimulus path for each existing stimulus
    % For the runs while this was made, there were 12 muscles in the sim, 12 dc stims, and 12 emg waveforms

    
    % Convert the emg signals to activation
    activation = emg2activation(inEMGdata{1});
    activation(activation<0) = 0;
    
    % Define the initial point, upper, and lower bounds
        % The upper bound is limited such that the maximum values of the activation do not exceed 20nA
        % ub is the same as "20./max(20*activation)"
    ub = 1./max(activation);
    ub(isinf(ub))=0;
    lb = zeros(size(ub));
    
    if nargin < 3 || isempty(initialPoint)
        initialPoint = ub./2;
    end

%     pattOpts = optimoptions('patternsearch','Display','iter','PlotFcn',{'psplotbestx','psplotbestf','psplotmeshsize'},...
%         'UseParallel',true,'MaxTime',600);
    pattOpts = optimoptions('patternsearch','UseParallel',true,'MaxTime',2*60,'Display','iter','MeshExpansionFactor',4,'InitialMeshSize',1);
    fun = @(scaleVec) objFun(simPath,scaleVec,activation,jointMotion);
    fun(initialPoint)
    [outScale,fVal] = patternsearch(fun,initialPoint,[],[],[],[],lb,ub,[],pattOpts);

    function [outVal,jointMotion] = objFun(simPath,scaleVec,activation,normalMotion)
        % Create a directory for this thread to work in
        [~,simName] = fileparts(simPath);
        [jobDir,jobID] = fileparts(tempname);
        clusterDir = [jobDir,'\',jobID];
        mkdir(clusterDir);
        % Copy an existing simulation file into that directory
        copyfile(simPath, clusterDir);
        jobSimPath = [clusterDir,'\',jobID,'.asim'];
        % Rename that .asim file using a temporary name
        movefile([clusterDir,'\',simName,'.asim'], [clusterDir,'\',jobID,'.asim']);

        % Import the new .asim file and modify the Joint Motion output file name
        inText = importdata(jobSimPath);
        txtInd = find(contains(inText,'OutputFilename>Joint'));
        inText{txtInd} = replaceBetween(inText{txtInd},'ion','.txt',['_',jobID]);

        % Generate the direct current files for each stimulus/muscle
        nsys = struct();
        nsys.proj_params.physicstimestep = .54;
        nsys.proj_params.simendtime = 10.01;
        for ii = 1:12
            %dataWave = scaleVec(ii)*[activation(:,ii);activation(:,ii);activation(:,ii)]*20;
            dataWave = scaleVec(ii)*activation(:,ii)*20;
            dcStimPath = [clusterDir,'\',jobID,'_',num2str(ii),'.txt'];
            generate_direct_current_file(nsys,dataWave,dcStimPath);
            stimDocInd = find(contains(inText,['<ID>stDC1-neur',num2str(ii),'-ID</ID>']))+11;
            inText{stimDocInd} = replaceBetween(inText{stimDocInd},'>','<',dcStimPath);
        end

        % Save the modifications we've made to the .asim file to a "mod" file
        jobSavePath = [jobSimPath(1:end-5),'_mod',jobSimPath(end-4:end)];
        fileID = fopen(jobSavePath,'w');
        fprintf(fileID,'%s\n',inText{:});
        fclose(fileID);

        % Run the mod .asim file using the system simulation execution
        jobSavePath = [jobSimPath(1:end-5),'_mod',jobSimPath(end-4:end)];
        sour_folder = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin';
        executable = [string([sour_folder,'\AnimatSimulator']),string(jobSavePath)];
        jsystem(executable);

        % Import the Joint Motion data from the output .txt file
        jointMotionPath = [clusterDir,'\','JointMotion_',jobID,'.txt'];
        ds = importdata(jointMotionPath);
        jointMotion = ds.data(:,3:5);
        
            temp(:,1) = jointMotion(:,find(contains(ds.colheaders,'Hip'))-2);
            temp(:,2) = jointMotion(:,find(contains(ds.colheaders,'Knee'))-2);
            temp(:,3) = jointMotion(:,find(contains(ds.colheaders,'Ankle'))-2);
            jointMotion = temp.*(180/pi)+[98.4373 102.226 116.2473];
            
        % Assign the output value of the objective function
        jointGrades = vecnorm(normalMotion(1:length(jointMotion),:)-jointMotion,2,1);
        %outVal = 10*jointGrades(1)+10*jointGrades(2)+jointGrades(3);
        outVal = sum(jointGrades.^2);
        
        % Saving this display code for convenience, can delete later
%         figure;subplot(3,1,1);
%             plot(ds.data(:,2),jointMotion(:,1),'LineWidth',2);hold on;plot(ds.data(1:length(jointMotion),2),normalMotion(1:length(jointMotion),1),'LineWidth',2);
%             xlim([0 10]); title('Hip');ylabel('Joint Angle (deg)');xlabel('Time (s)');
%             subplot(3,1,2)
%             plot(ds.data(:,2),jointMotion(:,2),'LineWidth',2);hold on;plot(ds.data(1:length(jointMotion),2),normalMotion(1:length(jointMotion),2),'LineWidth',2);
%             xlim([0 10]); title('Knee');ylabel('Joint Angle (deg)');xlabel('Time (s)');
%             subplot(3,1,3)
%             plot(ds.data(:,2),jointMotion(:,3),'LineWidth',2);hold on;plot(ds.data(1:length(jointMotion),2),normalMotion(1:length(jointMotion),3),'LineWidth',2);
%             xlim([0 10]); title('Ankle');ylabel('Joint Angle (deg)');xlabel('Time (s)');

        % Delete all files from the directory and the directory itself
        % This is time consuming but prevents massive storage issues from long runs (sometimes on the order of GBs)
        delete([clusterDir,'\*'])
        rmdir(clusterDir)
    end
end