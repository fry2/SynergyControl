%Sim Ani RMS finds the root mean square of the kinematics for both animal and
%simulation data and plots them both. Does only one front leg and one hind
%leg at a time.

function [AllAnimalMeanFront,AllAnimalMeanHind,AnimalStrideBackLegs] = Sim_Ani_RMS(to_plot, Stride, NewTime)
    close all

    %Define a time vector for interpolation of all steps to the same length
    Time = 0:.001:1;

    if to_plot
        %Plot everything on one graph with multiple subplots
        figure
        hold on
    end


    %% Hind legs
    %Cycle through all 3 joints on the hind legs.
    for k=1:3
        %% The Animal Section
        clear AnimalStride AnimalMean AnimalStd SimStride SimMean SimStd X Y
        % length(Stride)

        %Store the interpolated strides into a vector. Use 3 or 4 for hind left or
        %right leg.
        for i=1:length(Stride)
            AnimalStride(:,i) = interp1(NewTime{i,3},Stride{i,3}(:,k),Time);
        end

        %Find the mean value and standard deviation for all strides
        for j=1:size(AnimalStride,1)
            clear temp
            temp=AnimalStride(j,:);
            AnimalMean(j) = mean(temp(~isnan(temp)),2);
            AnimalStd(j) = std(temp(~isnan(temp)));
        end
        
        if to_plot
            %Plot the data with the standard deviation around it. The fliplr stuff is a
            %clever way of making the filled in box.
            X=[Time,fliplr(Time)];
            Y=[AnimalMean+AnimalStd,fliplr(AnimalMean-AnimalStd)];
            subplot(4,2,2*k+1)
            hold on
            fill(X,Y,[.5 .5 1])
            plot(Time,AnimalMean,'--b','Linewidth',2)
        end

        %Save the mean of the joint for later
        AllAnimalMeanHind(:,k) = AnimalMean';
        AnimalStrideBackLegs(:,:,k) = AnimalStride;
        %% The simulation section

    %     %Store the interpolated strides into a vector. Use 3 or 4 for hind left or
    %     %right leg
    %     for i=1:length(StrideS)
    %         if isempty(StrideS{i,4})==0
    %             SimStride(:,i) = interp1(NewTimeS{i,4},StrideS{i,4}(:,k),Time);
    %         else
    %             SimStride(:,i) = Time*NaN; 
    %         end
    %     end
    % 
    %     %Finds the mean value and standard deviation for all strides
    %     for j=1:size(SimStride,1)
    %         clear temp
    %         temp=SimStride(j,:);
    %         SimMean(j) = mean(temp(~isnan(temp)),2);
    %         SimStd(j) = std(temp(~isnan(temp)));
    % 
    %     end
    % 
    %     %Store the rms value for possible comparisons later
    %     HindRMS(k) = rms(AnimalMean-SimMean);
    % 
    %     %Plot the data with the standard deviation around it. The fliplr stuff is a
    %     %clever way of making the filled in box.
    %     X=[Time,fliplr(Time)];
    %     Y=[SimMean+SimStd,fliplr(SimMean-SimStd)];
    %     subplot(4,2,[2*k+1])
    %     hold on
    %     fill(X,Y,[.5 1 1])
    %     plot(Time,SimMean,'k','Linewidth',2)

    end

    %% Front legs, the exact same as hind legs, see comments above

    for k=1:4

        clear AnimalStride AnimalMean AnimalStd SimStride SimMean SimStd  X Y

        for i=1:length(Stride)
            if isempty(Stride{i,2})==0
                AnimalStride(:,i) = interp1(NewTime{i,2},Stride{i,2}(:,k),Time);
            else
               AnimalStride(:,i) = Time*NaN; 
            end
        end

        for j=1:size(AnimalStride,1)
            clear temp
            temp=AnimalStride(j,:);
            AnimalMean(j) = mean(temp(~isnan(temp)),2);
            AnimalStd(j) = std(temp(~isnan(temp)));
        end

        if to_plot
            % figure
            % hold on
            X=[Time,fliplr(Time)];
            Y=[AnimalMean+AnimalStd,fliplr(AnimalMean-AnimalStd)];
            subplot(4,2,2*k)
            hold on
            fill(X,Y,[.5 .5 1])
            plot(Time,AnimalMean,'--b','Linewidth',2)
        end

        AllAnimalMeanFront(:,k) = AnimalMean';

    %     for i=1:length(StrideS)
    %         if isempty(StrideS{i,2})==0
    %             SimStride(:,i) = interp1(NewTimeS{i,2},StrideS{i,2}(:,k),Time)';
    %         else
    %             SimStride(:,i) = Time*NaN; 
    %         end
    %     end
    % 
    %     for j=1:size(SimStride,1)
    %         clear temp
    %         temp=SimStride(j,:);
    %         SimMean(j) = mean(temp(~isnan(temp)),2);
    %         SimStd(j) = std(temp(~isnan(temp)));
    %     end
    % 
    %     % FrontRMS(k) = rms(AnimalMean-SimMean);
    % 
    %     X=[Time,fliplr(Time)];
    %     Y=[SimMean+SimStd,fliplr(SimMean-SimStd)];
    %     subplot(4,2,2*k)
    %     hold on
    %     fill(X,Y,[.5 1 1])
    %     plot(Time,SimMean,'k','Linewidth',2)


    end
end

