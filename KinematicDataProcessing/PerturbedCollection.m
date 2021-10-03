%Perturbed Collection loads in simulation data from the declared folder
%inside and for the number of trials desired. It then does a similar
%normalization process to that done with the animal data. Then it finds the
%mean stride and standard deviation of all the steps in the desired time
%period for all the trials.
function [SFootContact, SFrontLeft, SFrontRight, SBackLeft, SBackRight, STime] = PerturbedCollection(to_plot)
    close all
    clear AllSCoordinationMean AllSCoordinationSTD tempmean tempstd StepRecov StepRecovAvg idcount Steps
    %%
    idcount=1;

    %the id declairs what trials to use, can be a vector of ever increasing
    %values
    for id=2
        if id==1
            DataToGet1 = 'DataTool_1.txt';
            DataToGet2 = 'DataTool_2.txt';
        else
            DataToGet1 = ['DataTool_1 (',sprintf(num2str(id)),').txt'];
            DataToGet2 = ['DataTool_2 (',sprintf(num2str(id)),').txt'];
        end


    %This is the folder to look in for the trials
    folder = [fileparts(mfilename('fullpath')),'\UnPerturbed\NewInter\'];

    Path1 = [folder,DataToGet1];
    Path2 = [folder,DataToGet2];

    %Process Simulation Kinematics loads in the simulation data from the above
    %folder and organize it for further fun
    [SFootContact, SFrontLeft, SFrontRight, SBackLeft, SBackRight, STime] = Process_Simulation_Kinematics(to_plot, Path1,Path2);

    %Simulation Kinematics plots the foot contact vs times. Inside is a start
    %and end time for cropping the data for a closer look plot. This function
    %does nothing for further functions.
    Simulation_Kinematics(to_plot, SFootContact, STime);

    %Sim Coodination finds the timing as a percent of stride for foot
    %touchdown. Inside has a beginpoint and stoppoint for telling the system
    %where to perform the analysis within the trial. It does similar step
    %timing and normalization as is done inside Animal Coordination. The major
    %difference is that there is included a buffer time which does not allow
    %for steps to be shorter than .14 seconds. Does not use kinematics.
    [SCoordinationMean,SCoordinationSTD] = Sim_Coordination(SFootContact);
    % keyboard

    %Convolve Sim is used to find the time it takes to recover from the
    %perturbation (if there is one). It uses a range of time before the
    %perturbation, then conlves it with after the perturbation to find when
    %step timing is similar again.
    [StepRecov(idcount),temp] = ConvolveSim(to_plot, STime,SFootContact);

    %Store the coordination mean and standard deviation of the individual trial
    %away.
    AllSCoordinationMean(:,:,idcount) = SCoordinationMean;
    AllSCoordinationSTD(:,:,idcount) = SCoordinationSTD;

    idcount=idcount+1;

    end
    %%
    %Find Mean of all the trials' means and standard deviation of all the
    %trials' standard deviation
    for id=1:4
        for j=1:4
            clear temp
            tempmean = nonzeros(AllSCoordinationMean(id,j,:));
            tempstd = nonzeros(AllSCoordinationSTD(id,j,:));
                if id+j==5 || id+j==5 || id+j==7
                    tempmean=tempmean-[tempmean>.5];
                end
            CoordinationMean(id,j) = mean(tempmean);
            CoordinationSTD(id,j) = sqrt(sum(tempstd.^2/length(tempstd)));
        end
    end

    %Show us the results
%     CoordinationMean
%     CoordinationSTD
    StepRecov=floor(StepRecov);
    finiteStep = StepRecov(isfinite(StepRecov));
    StepRecovAvg = mean(StepRecov(isfinite(StepRecov)));
    Steps = [max(finiteStep),min(finiteStep),StepRecovAvg,sum(isinf(StepRecov))];
    save([folder,'Mean STD2 and steps'],'CoordinationMean','CoordinationSTD','Steps')
end