%Animal Stride normalizing uses the animal kinematic and ground contact
%data to determine when each step starts and when each step ends. It then
%uses these beginings and endings to set stance at 0-50% of the stride and
%swing at 51-100%. It then plots out all the joints.
function [Stride, StrideContact, NewTime] = Animal_Stride_Normalizing(to_plot, ATime, AFootContact, AFrontLeft, AFrontRight, ABackLeft, ABackRight)
    %Clear variables which are rewritten in this subpart of the program
    clear Stride StrideTime StrideContact NewTime Timelength Timing Allstrides
    close all

    %Cycle through all legs
    for k=1:4
        as = 1;
        %Cycle through trials (leaving out fastest and slowest)
        for i=1:8

            clear TData TFootContact start delaystart delayend endstride middle

            %Initialize Variables with data from the correct leg
            if k==1
                TAData = AFrontLeft{i};
            elseif k==2
                TAData = AFrontRight{i};
            elseif k==3
                TAData = ABackLeft{i};
            elseif k==4
                TAData = ABackRight{i};
            end
            TAFootContact = AFootContact{i}(:,k);
            start=0;
            delaystart=0;
            delayend=1;
            endstride=0;

            %Check if foot is starting on ground. if it is don't make that start of
            %stride
            if TAFootContact(1)<0
                delaystart=1;
            end

            %Loop to find beginning and end of stride (start with foot contact)
            %only bothers finding one stride per leg
            for j=2:length(TAFootContact)
                %Special case (maybe) needed for trial number 9
            %         if i == 9
            %             if TAData(j)>360
            %                 TAData(j) = TAData(j)-360
            %             end
            %         end
                %Only choose a start if foot doesn't start on ground
                    if delaystart==1
                        if isnan(TAFootContact(j))==1
                            delaystart=0;
                        end
                    else
                        if TAFootContact(j)<0 && start==0
                            start=j;
                        end
                    end

                %Choose an end after start has been found and foot goes off ground
                    if delayend==1 && start>0
                        if isnan(TAFootContact(j))==1
                            middle=j;
                            delayend=0;
                        end
                    else
                        if TAFootContact(j)<0 && start~=0 && endstride==0
                            endstride=j-1;
                        end
                    end

            end

            %Start = foot touches the ground
            %Middle = foot just lifted from the ground
            %Endstride = foot touches ground again

            %After going through the entire vector, if we found a full stride,
            %then save it. The animal trials are so short we will only be able
            %to find one full stride and don't worry about repeats.
            if endstride~=0
                Stride{i,k} = TAData(start:endstride,:);
                StrideTimeStance{i,k} = ATime{i}(start:middle);
                StrideTimeSwing{i,k} = ATime{i}(middle+1:endstride);
                StrideContact{i,k} = TAFootContact(start:endstride);
            end
        end

        %Normalize time between 0 and 1
        for i=1:length(Stride)
            if isempty(Stride{i,k})==0        
                TimelengthStance = StrideTimeStance{i,k}(end)-StrideTimeStance{i,k}(1);
                TimelengthSwing = StrideTimeSwing{i,k}(end)-StrideTimeSwing{i,k}(1);
                TimingStance = StrideTimeStance{i,k};
                TimingSwing = StrideTimeSwing{i,k};
                NewTime{i,k} = [(TimingStance - TimingStance(1))/(2*TimelengthStance+.001),.5+(TimingSwing - TimingSwing(1))/(2*TimelengthSwing)];
            end
        end
        
        if to_plot
            %% Plot all of the Joints with Ground Contact underneath
            %Hind leg
            if k==3 || k==4
                figure(1)
                hold on
                figure(2)
                hold on
                figure(3)
                hold on

                for i=1:length(Stride)
                    if isempty(Stride{i,k})==0
                        figure(1)
                        plot(NewTime{i,k},Stride{i,k}(:,1),'b','LineWidth',2)
                        figure(2)
                        plot(NewTime{i,k},Stride{i,k}(:,2),'g','LineWidth',2)
                        figure(3)
                        plot(NewTime{i,k},Stride{i,k}(:,3),'r','LineWidth',2)

                        %Store the stride and time away
                        Allstrides{as,k} = Stride{i,k};
                        GoodTimes{as,k} = NewTime{i,k};
                        as = as+1;
                    end
                end
            end

            %Front Leg
            if k==1 || k==2
                figure(4)
                hold on
                figure(5)
                hold on
                figure(6)
                hold on
                figure(7)
                hold on

                %Plot the stride if it exists
                for i=1:length(Stride)
                    if isempty(Stride{i,k})==0
                        figure(4)
                        plot(NewTime{i,k},Stride{i,k}(:,1),'b','LineWidth',2)
                        figure(5)
                        plot(NewTime{i,k},Stride{i,k}(:,2),'m','LineWidth',2)
                        figure(6)
                        plot(NewTime{i,k},Stride{i,k}(:,3),'c','LineWidth',2)
                        figure(7)
                        plot(NewTime{i,k},Stride{i,k}(:,4),'k','LineWidth',2)

                        %Store the stride and time away
                        Allstrides{as,k} = Stride{i,k};
                        GoodTimes{as,k} = NewTime{i,k};
                        as = as+1;
                    end
                end
            end
        end
    end

    if to_plot
        %Add titles and labels to the figures
        legendvals={'Hip';'Knee';'Ankle';'Scapula';'Shoulder';'Elbow';'Wrist'};
        for i=1:7
            figure(i)
            title(legendvals{i});
            xlabel('Stride Fraction');
            ylabel('Angle (Deg)');
        end
    end

    %Find how long the swings occurs
    swings = reshape(StrideTimeSwing,1,size(StrideTimeSwing,1)*4);
    a = 1;
    for i=1:length(swings)
        if ~isempty(swings{i})
            len(a) = length(swings{i});
            a = a+1;
        end
    end

    %Multiply length of vector by timestep
    mean_time_of_swing = mean(len)*.002;
end