%Simulation Kinematics plots the foot contact vs times. Inside is a start
%and end time for cropping the data for further analysis.
function Simulation_Kinematics(to_plot, SFootContact, STime)
    clear SFootContactPlot
    if to_plot
        hfig=figure;
        % hold on

        for i=1:4
            SFootContactPlot(:,i)=SFootContact(:,i)-.03*i;
        end

        % lengthfill = ones(1,(32500-30000)/2)';


        % plot([3 3],[.8501 1],':k')
        % plot([3.25 3.25],[.8501 1],':k')
        % plot(STime(20000/2+1:45000/2),SFootContactPlot(20000/2+1:45000/2,:),'Linewidth',15)
        % set(hfig, 'Position', [500 500 400 150])

        % lengthfill = ones(1,10000/2)';

        %What is the start and end time you want to plot in seconds
        Cstart=1;
        Cend=10;

        %Makes the figure with the correct limits
        axes1 = axes('Parent',hfig);
        set(hfig, 'Position', [50 300 500 150])
        xlim(axes1,[0 Cend-Cstart]);
        hold(axes1,'all');

        %Plot the foot contacts
        CTime = STime(Cstart*10000/2+1:Cend*10000/2)-Cstart;
        CFootContactPlot=SFootContactPlot(Cstart*10000/2+1:Cend*10000/2,:);
        % plot([3 3],[.8501 1],':k')
        % plot([3.25 3.25],[.8501 1],':k')
        plot(CTime,CFootContactPlot,'Linewidth',15)

        figure(hfig)
    end
end
