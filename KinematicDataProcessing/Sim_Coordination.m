%Sim Coodination finds the timing as a percent of stride for foot
%touchdown. Inside has a beginpoint and stoppoint for telling the system
%where to perform the analysis within the trial. It does similar step
%timing and normalization as is done inside Animal Coordination. The major
%difference is that there is included a buffer time which does not allow
%for steps to be shorter than .14 seconds.

function [SCoordinationMean,SCoordinationSTD] = Sim_Coordination(SFootContact)

    clear start AllPercents lengthstarts CurrentTrial lengthtrial numtrials AllStarts AllEnds CompactPercents ACoordinationSTD CoordinationMean

    % numtrials = size(AFootContact,2);

    %What is the begin and end points to be analyzed? 30000/2+1 is 3 seconds in
    %as the time step is .2 ms
    beginpoint=30000/2+1;
    stoppoint=90000/2;

    %Cycle through each foot
    for i=1:4    
        delay=0;
        footdown=0;
        nexti=1;
        CurrentFoot = SFootContact(beginpoint:stoppoint,i);

        %Check if foot is starting on ground. if it is don't make that start of
        %stride   
        if CurrentFoot(1)==1
            delay=1;
            footdown=1;
        end

        %Loop to find beginning and end of stride (start with foot contact)
        for j=151:length(CurrentFoot)
            %If foot started on the ground, figure out when it gets off the ground
            if delay==1
                if isnan(CurrentFoot(j))&& isnan(CurrentFoot(j-150));
                    delay=0;
                    footdown=0;
                end
            end

            if delay==0 %If the foot has been off the ground
                %Find liftoff points as long as it is greater than .14 seconds from
                %touchdown
                if footdown==1 && (j-start(nexti,i))>1400
                    if isnan(CurrentFoot(j));
                        footdown=0;
                        endstride(nexti,i)=j;
                        nexti=nexti+1;
                    end
                elseif footdown==0
                    %Keep track of everytime the foot hits the ground again
                    if CurrentFoot(j)==1
                        start(nexti,i)=j;
                        footdown=1;
                    end
                end
            end
        end
    end

    %Find the length of each stride in number of datapoints
    for i=2:size(start,1)
        if start(i,:)~=0
            lengthstarts(i-1,:) = start(i,:)-start(i-1,:);
        end
    end

    %Process touchdown and liftoffs to get a percent stride
    for i=1:4
        clear percent CurrentStarts
        CurrentStarts = start;

        %Compare other legs to this leg
        for j=1:4

            %Compare for each full stride
            for k=1:size(CurrentStarts,1)-1

                %Only check if there is a full stride
                if CurrentStarts(k,j)~=0 && CurrentStarts(k+1,j)~=0
                    stridelength = CurrentStarts(k+1,j)-CurrentStarts(k,j);

                    %Compare all other legs to leg j
                    for l=1:4

                        if l~=j
                            %Compare all possible touchdowns
                            for m=1:size(CurrentStarts,1)-1

                                %Only check if there is a touchdown
                                TouchdownTime = CurrentStarts(m,l)-CurrentStarts(k,j);
                                realpercent = TouchdownTime/stridelength;
                                if CurrentStarts(m,l)~=0 && TouchdownTime>0 && realpercent<1
                                    percent(k,m,l,j) = realpercent;
                                end
                            end
                        end
                    end
                end
            end
        end

        AllPercents = percent;
    end

    %Gather the data more compactly
    %Initial touchdown leg
    for i=1:4
        %Following Leg
        for j=1:4  
            nexti=1;
            %Each Trial
            for l=1:1
                clear CurrentPercent
                CurrentPercent = AllPercents;
                %Each stride
                if l~=7 && (j~=4 || i~=4)
                    for k=1:size(CurrentPercent,1)
                        %Each Compare
                        for m=1:size(CurrentPercent,2)
                            if CurrentPercent(k,m,i,j)~=0
                                CompactPercents(i,j,nexti)=CurrentPercent(k,m,i,j);
                                nexti=nexti+1;
                            end
                        end
                    end
                end
            end
        end
    end

    %Find the mean and standard deviation of the percentages for leg timing for
    %certain values, if it is greater than 50%, make it a negative value
    %instead (aka -.1% instead of .9%)
    for i=1:4
        for j=1:4
            clear temp
            temp = nonzeros(CompactPercents(i,j,:));
                if i+j==5
                    temp=temp-[temp>.5];
                end
            SCoordinationMean(i,j) = mean(temp);
            SCoordinationSTD(i,j) = std(temp);
        end
    end
end
