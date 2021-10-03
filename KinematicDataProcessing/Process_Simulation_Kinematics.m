%Load data from Animatlab Simulation and organize it

%clear SImportKin SData Stime SFrontLeft SFrontRight SBackLeft SImportContact ScontactData SFootContact
function [SFootContact, SFrontLeft, SFrontRight, SBackLeft, SBackRight, STime] = Process_Simulation_Kinematics(to_plot, Path1,Path2)
% Path1 =  'C:\Users\Hunt\Desktop\PhD\AnimatLab2\Rat 14-7-11 - Unsupported hind leg walking\Angle Comp.txt';
% Path2 =  'C:\Users\Hunt\Desktop\PhD\AnimatLab2\Rat 14-7-11 - Unsupported hind leg walking\Ft Pres.txt'';

%Load in data
SImportKin = importdata(Path1);

%Store data in an array
SData = SImportKin.data;

%Timestamp in first column
STime = SData(:,1);

%Store data into proper places and organize so that Hip Knee Ankle and
%Scapula Shoulder Elbow Wrist are the order, Scapula and Wrist are negative
%Left front joint angle data Original Order Elbow, Scapula, Shoulder, Wrist
SFrontLeft = [SData(:,9)+pi/4,SData(:,10)+3*pi/4,SData(:,8)+pi/2,SData(:,11)+pi]*180/pi;

%Right front joint angle data Original Order: Wrist Shoulder Scapula Elbow
SFrontRight = [SData(:,13)+pi/4,SData(:,14)+3*pi/4,SData(:,12)+pi/2,SData(:,15)+pi]*180/pi;

%Left back joint angle data Original Order: Knee Ankle Hip
SBackLeft = [SData(:,4)+pi/2,SData(:,3)+pi,SData(:,2)+pi/2]*180/pi;

%Right back joint angle data Original Order: Knee Ankle Hip
SBackRight = [SData(:,5)+pi/2,SData(:,7)+pi,SData(:,6)+pi/2]*180/pi;

%Import contact data
SImportContact = importdata(Path2);

SContactData = SImportContact.data;

%Remove timestamp and reorganize to be FL, FR, BL, BR
SFootContact(:,1:2) = SContactData(:,4:5);
SFootContact(:,3:4) = SContactData(:,2:3);

%Normalize contact data to read NaNs when not in contact, and 1 when in contact, similar to the Animal data
SFootContact(:,1:2) = SFootContact(:,1:2)>-.0602;
SFootContact(:,3:4) = SFootContact(:,3:4)>-.0599;

if to_plot
    figure
    plot(SFootContact(:,3),'b')
end

%Pad foot contact to make a little longer to make up for setting values
%higher
padamount = 50;
SFootContact(1:end-padamount,3:4) = SFootContact(1+padamount:end,3:4);

for i=length(SFootContact):-1:1+padamount*2
    if SFootContact(i-2*padamount,3) == 1
        SFootContact(i,3) = 1;
    end
    if SFootContact(i-padamount*2,4) == 1
        SFootContact(i,4) = 1;
    end
end

if to_plot
    hold on
    plot(SFootContact(:,3)+.01,'r')
end

SFootContact(SFootContact==0)=NaN;
end