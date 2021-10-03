% Runs through a selection of lambda values for the ACLS NNMF algorithm from Langville 2014
% Select a number of synergies
% Output meanMat has cols=lamda W and rows=lambda H
% If output on NMFdecomposition has changed, setup such that: [~,~,~,~,differ] = NMFdecomposition(k,forces,0,lamW,lamH);
% Notably, set the above function to take lamW and LamH as inputs
% Where differ is an array abs(inputForces-recombinedForces)

close all

lamHvals = linspace(.001,.1,20);
lamWvals = lamHvals;
meanMat = zeros(length(lamHvals),length(lamHvals));

tstart = tic;
for i = 1:length(lamHvals)
    lamH = lamHvals(i);
    rowHold = zeros(1,length(20));
    for j = 1:length(lamWvals)
        lamW = lamWvals(j);
        [~,~,~,~,differ] = NMFdecomposition(5,forces,lamW,lamH,0);
        meanDiff = mean(mean(differ));
        rowHold(j) = meanDiff;
    end
    meanMat(i,:) = rowHold;
end
telapsed = toc(tstart);
disp([num2str(telapsed),' s'])

figure
subplot(2,1,1)
plot(lamWvals,mean(meanMat))
ylabel('lamW')
subplot(2,1,2)
plot(lamHvals,mean(meanMat,2))
ylabel('lamH')