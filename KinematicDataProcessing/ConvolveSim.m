%Convolve Sim is used to find the time it takes to recover from the
%perturbation (if there is one). It uses a range of time before the
%perturbation, then conlves it with after the perturbation to find when
%step timing is similar again.

function [StepRecoveryat,SimCorr] = ConvolveSim(to_plot, STime2,SFootContact2)
    %Initiate the ranges of before and after perturbation to compare. After
    %should be longer than before
    range1start = 25000/2+1;
    range1end = 30000/2;
    range1=range1end-range1start;
    range2start = 32500/2;
    range2end = 50000/2;

    %Turn the contact timing into a double (from boolean?)
    SFootContact2 = double(~isnan(SFootContact2));

    %Grab the vectors of before and after
    BeforeP = SFootContact2(range1start:range1end,:);
    AfterP = SFootContact2(range2start:range2end,:);

    %Cross correlated the before time to amke sure it had actually found
    %coordination
    FullCorr = xcorr2(BeforeP);

    %Cross coorelate the before and after times
    SimCorr = xcorr2(BeforeP,AfterP);
    RelCorr = fliplr(SimCorr(range1+1:end,4)');

    %Find the points of highest correlation, making sure at least a specific
    %amount of time between peaks has passed, that is basically a full stride.
    %Value found through testing
    [peaks,Steps] = findpeaks(RelCorr,'Minpeakdistance',1500);

    %Determine a minimum correlation necessary for us to have declared
    %coordination is regained. Set at 80% for initial work
    numSteps = length(Steps);
    RelCorr = RelCorr - max(FullCorr(:,4))*.8;

    %Figure out when the 80% is first achieved and determine how many steps it
    %took after the perturbation 
    jd=1;
    TheRecovery=-1;
    while TheRecovery==-1
        Stepsjd=Steps(jd);
        Reljd=RelCorr(Steps(jd));

        if RelCorr(Steps(jd))>0
            if RelCorr(Steps(jd-1))>0
                TheRecovery=jd-1;
            end
        end

        jd=jd+1;

        %If coordination is not regained, fill in Inf
        if jd>numSteps && TheRecovery==-1
            TheRecovery=Inf;
        end
    end

    %Store the number of steps
    StepRecoveryat = TheRecovery;

    if to_plot
        %Plot the correlation for visual analysis
        figure
        plot(STime2(range2start:range2end),RelCorr)
    end
end

