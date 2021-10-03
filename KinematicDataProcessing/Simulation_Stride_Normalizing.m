%Normalize the data into a stride percentage, and break out all steps
%individually

function Simulation_Stride_Normalizing(to_plot, SFrontLeft, SFrontRight, SBackLeft, SBackRight, SFootContact, STime)
    clear StrideS StrideSTime StrideSContact NewTimeS Timelength Timing nexti StrideSStance StrideSSwing StrideSTimeStance StrideSTimeSwing TimelengthStance TimelengthSwing TimingStance TimingSwing
    % close all

    %Cycle through all legs
    for k=1:4

        clear TSData TSFootContact start delaystart delayend endstride beginpoint stoppoint

        beginpoint=20000/2;
        stoppoint=100000/2;

        %Initialize Variables with data
        if k==1
            TSData = SFrontLeft(beginpoint:stoppoint,:);
        elseif k==2
            TSData = SFrontRight(beginpoint:stoppoint,:);
        elseif k==3
            TSData = SBackLeft(beginpoint:stoppoint,:);
        elseif k==4
            TSData = SBackRight(beginpoint:stoppoint,:);
        end    
        %Initialize Variables with data
        TSFootContact = SFootContact(beginpoint:stoppoint,k);
        start=0;
        delaystart=0;
        nexti=1;
        endstride=0;

        %Check if foot is starting on ground. if it is don't make that start of
        %stride
        if TSFootContact(1)<0
            delaystart=1;
        end

        %Loop to find beginning and end of stride (start with foot contact)
        for j=2:length(TSFootContact)

            %Only choose a start if foot doesn't start on ground
                if delaystart==1
                    if isnan(TSFootContact(j))==1
                        delaystart=0;
                        endstride(nexti)=j;
                         nexti=nexti+1;
                    end
                else
                    %Keep track of everytime the foot hits the ground again
                    if TSFootContact(j)==1
                        start(nexti)=j;
                        delaystart=1;
                    end
                end       
        end

        %For each time we have a foot touchdown (minus 1) store the stride
        %components into the proper vectors and store the time vectors as well
        %for use in normalizing later (by adjusting the time vector)
        for l=3:length(start)-1
            if start(l)-start(l-1)>10 %Make sure the starts didn't happen too close together
                StrideSStance{l-1,k} = TSData(start(l-1):endstride(l-1),:);
                StrideSSwing{l-1,k} = TSData(endstride(l-1)+1:start(l),:);
                StrideSTimeStance{l-1,k} = STime(start(l-1):endstride(l-1));
                StrideSTimeSwing{l-1,k} = STime(endstride(l-1)+1:start(l));
                StrideS{l-2,k} = [StrideSStance{l-1,k};StrideSSwing{l-1,k}];
                StrideSContact{l-2,k} = TSFootContact(start(l-1):start(l));
            end
        end

        %Normalize time between 0 and 1 with stance in the first half and swing
        %in the second half
        for l=2:length(StrideS)+1
            if isempty(StrideS{l-1,k})==0        
                TimelengthStance = StrideSTimeStance{l,k}(end)-StrideSTimeStance{l,k}(1);
                TimingStance = StrideSTimeStance{l,k};
                TimelengthSwing = StrideSTimeSwing{l,k}(end)-StrideSTimeSwing{l,k}(1);
                TimingSwing = StrideSTimeSwing{l,k};
                NewTimeS{l-1,k} = [(TimingStance - TimingStance(1))/(2*TimelengthStance+.001);.5+(TimingSwing - TimingSwing(1))/(2*TimelengthSwing)];
            end
        end

        if to_plot
            %Plot all of the Joints with Ground Contact underneath
            if k==1 || k==3
            figure
            hold on
            end
            for i=1:length(StrideS)
                plot(NewTimeS{i,k},StrideS{i,k},NewTimeS{i,k},StrideSContact{i,k}+2*i-100)
            end
        end

    end
end