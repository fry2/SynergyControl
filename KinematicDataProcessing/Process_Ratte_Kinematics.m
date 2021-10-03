%Process Ratte Kinematics loads in the data from Ratte 1 Limb
%kinematics.xlsx
%which we recieved from Manuela Schmidt giving kinematic traces and ground
%contact timing of 10 rat trials. It then organizes it for further
%processing and comparison with other functions.

function [ATime, AFootContact, AFrontLeft, AFrontRight, ABackLeft, ABackRight] = Process_Ratte_Kinematics

    clear all
    close all

    %Loop to read in data from all sheets
    for i=1:10

    %Clear out reused variables
    clear Num Txt Raw Data

    %Read in data
    [Num,Txt,Raw] = xlsread('G:\My Drive\Rat\SynergyControl\Data\Ratte 1 Limb kinematics.xlsx',i);

    %Velocity from cell B2
    Avelocity(i) = Num(1,2);

    %Eliminate top 4 rows which do not contain any more data
    AData = Num(5:end,:);

    %Time in first column
    ATime{i} = AData(:,1)';

    %Time normalized to start at 0
    ATime{i} = ATime{i}-ATime{i}(1);

    %Contact in columns C-F. Given as random negative number when ground in
    %contact (useful for plotting) - Front Left, Front Right, Back Left, Back
    %Right
    AFootContact{i} = AData(:,3:6);

    %Left front joint angle data in columns G (Scap angle) K-M (Shoulder,
    %Elbow, Wrist)
    AFrontLeft{i} = [AData(:,7),AData(:,11:13)];

    %Right front joint angle data in columns N (Scap angle) R-T (Shoulder,
    %Elbow, Wrist)
    AFrontRight{i} = [AData(:,14),AData(:,18:20)];

    %Left back joint angle data in columns Y-AA (Hip, Knee, Ankle)
    ABackLeft{i} = AData(:,25: 27);

    %Right back joint angle data in columns AF-AH (Hip, Knee, Ankle)
    ABackRight{i} = AData(:,32:34);

    end
    clear AData i
    % 
    % Extra code for plotting things when loading it in to check the data
    %
    % figure
    % hold on
    % for i=1:10
    %     figure
    %     plot(ATime{i},AFrontLeft{i}(:,1),ATime{i},AFrontLeft{i}(:,2),ATime{i},AFrontLeft{i}(:,3),ATime{i},AFrontLeft{i}(:,4),ATime{i},AFootContact{i}(:,1))
    % end
    % 
    % figure
    % hold on
    % for i=1:10
    %     plot(ATime{i},AFootContact{i}(:,1),ATime{i},AFootContact{i}(:,2),ATime{i},AFootContact{i}(:,3),ATime{i},AFootContact{i}(:,4));
    % end

    %Plot foot contact data
    % figure
    % plot(Time,FootContact(:,1),Time,FootContact(:,2),Time,FootContact(:,3),Time,FootContact(:,4));
    % legend('Front Left','Front Right','Back Left','Back Right')
    % title('Foot contacts')
    % 
    % %Plot Scapula Angle with Left Front ground contact
    % figure
    % plot(Time,FrontLeft(:,1),Time,FootContact(:,1))
    % title('Left Scapula Angle')
    % 
    % %Plot Left Shoulder Angle with Left Front ground contact
    % figure
    % plot(Time,FrontLeft(:,2),Time,FootContact(:,1))
    % title('Left Shoulder Angle')
end