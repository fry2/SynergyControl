%Animal Coordination finds the mean and standard deviation for all the
%animal data for timing of footfall between different legs.

function [ACoordinationMean, ACoordinationSTD] = Animal_Coordination(AFootContact)
    clear AllPercents CurrentTrial lengthtrial numtrials AllStarts AllEnds CompactPercents ACoordinationSTD CoordinationMean

    numtrials = size(AFootContact,2);


    %First find all the foot touchdowns and liftoffs
    %Loop through all trials
    for k=1:numtrials

        clear start endstride

        start = ones(1,4)*NaN;
        CurrentTrial = AFootContact{k};
        lengthtrial = size(CurrentTrial,1);

        %Loop through all legs
        for i=1:4

            clear nexti delaystart 

            %% Find stride of leg

            %Initialize recording space
            nexti=1;

            %Determine if foot starts on ground
            if CurrentTrial(1,i)<0
                delaystart=1;
            else
                delaystart=0;
            end

            for j=2:lengthtrial
            %Only choose a start if foot doesn't start on ground
                if delaystart==1
                    if isnan(CurrentTrial(j,i))==1
                        delaystart=0;
                        %Keep track of liftoff
                        endstride(nexti,i)=j;
                        nexti=nexti+1;
                    end
                else
                    %Keep track of everytime the foot hits the ground again
                    if CurrentTrial(j,i)<0
                        start(nexti,i)=j;
                        delaystart=1;
                    end
                end       
            end


            AllStarts{k} = start;
            AllEnds{k} = endstride;

        end
    end

    %Process touchdown and liftoffs to get a percent stride
    for i=1:numtrials

        clear percent CurrentStarts
        CurrentStarts = AllStarts{i};

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

        %Store away the percent for each trial
        AllPercents{i} = percent;

    end

    %Gather the data more compactly
    %Initial touchdown leg
    for i=1:4
        %Following Leg
        for j=1:4  
            nexti=1;
            %Each Trial
            for l=1:numtrials
                clear CurrentPercent
                CurrentPercent = AllPercents{l};
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

    %Find the mean and standard deviation of the percents. For specific
    %positions, if the value is above .5, make it a negative % instead (aka
    %-.1% instead of .9%
    for i=1:4
        for j=1:4
            clear temp
            temp = nonzeros(CompactPercents(i,j,:));
            if i+j==5
                    temp=temp-[temp>.5];
                end
            ACoordinationMean(i,j) = mean(temp);
            ACoordinationSTD(i,j) = std(temp);
        end
    end
end