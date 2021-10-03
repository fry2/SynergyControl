function [Stride, StrideContact, NewTime, BackMean,BackRaw] = CompareAnimalAndSim(to_plot)
    %This file runs several functions\sections of code in order to load in and
    %process the animal data, simulation data, and make comparisons of the two
    %data sets.

    %Process Ratte Kinematics loads in the data from Ratte 1 Limb kinematics
    %which we recieved from Manuela Schmidt giving kinematic traces and ground
    %contact timing of 10 rat trials. It then organizes it for further
    %processing and comparison with other functions.
    [ATime, AFootContact, AFrontLeft, AFrontRight, ABackLeft, ABackRight] = Process_Ratte_Kinematics;

    %Animal Stride normalizing uses the animal kinematic and ground contact
    %data to determine when each step starts and when each step ends. It then
    %uses these beginings and endings to set stance at 0-50% of the stride and
    %swing at 50-100%. It then plots out all the joints.
    [Stride, StrideContact, NewTime] = Animal_Stride_Normalizing(to_plot, ATime, AFootContact,...
        AFrontLeft, AFrontRight, ABackLeft, ABackRight);

    %Perturbed Collection loads in simulation data from the declared folder
    %inside and for the number of trials desired. Then it finds the
    %mean stride timing and standard deviation of all the steps in the desired time
    %period for all the trials.
    [SFootContact, SFrontLeft, SFrontRight, SBackLeft, SBackRight, STime] = PerturbedCollection(to_plot);

    %Animal Coordination finds the mean and standard deviation for all the
    %animal data for timing of footfall between different legs.
    [ACoordinationMean, ACoordinationSTD] = Animal_Coordination(AFootContact);

    %Simulation Kinematics plots the foot contact vs times. Inside is a start
    %and end time for cropping the data for further analysis.
    Simulation_Kinematics(to_plot, SFootContact, STime);

    %Simulation Stide Normalizing is similar to Animal Stride Normalizing,
    %where each step is normalized such that stance is 0-50% and swing is
    %50-100%
    Simulation_Stride_Normalizing(to_plot, SFrontLeft, SFrontRight, SBackLeft, SBackRight, SFootContact, STime);

    %Sim Ani RMS finds the root mean square of the kinematics for both animal and
    %simulation data and plots them both. Does only one front leg and one hind
    %leg at a time.
    [AllAnimalMeanFront,BackMean,BackRaw] = Sim_Ani_RMS(to_plot, Stride, NewTime);
end