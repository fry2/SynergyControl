% Create correlation figures and plotWH

%[vafscores,recompiled,W,H,differ] = NMFdecomposition(6,forces,.01,.01,1);
%[outk,vafscores,recompiled,W,H] = NMFsyncounter(forces);

%% Check all open figures to see if planned figures are already open. Close all currently existing planned figures.
figHandles = get(groot, 'Children');
if ~isempty(figHandles)
    priorFigures = contains({figHandles(:).Name},{'wFig','hFig'});
    close(figHandles(priorFigures))
end

% numK = size(H,1);
% hMat = zeros(numK,numK);
% wMat = zeros(numK,numK);
% for dd = 1:numK
%     for ee = 1:numK
%         hMat(dd,ee) = corr(H(dd,:)',H(ee,:)','Type','Pearson');
%         wMat(dd,ee) = corr(W(:,dd),W(:,ee),'Type','Pearson');
%     end
% end

hMat = corr(H',H');
wMat = corr(W,W);

%% Create corrFigs
wFig = createCorrFig(wMat);
wFig.Children(2).Title.String = 'W Correlation';
set(wFig,'name','wFig','Position',[-1920,129,944,988])
hFig = createCorrFig(hMat);
hFig.Children(2).Title.String = 'H Correlation';
set(hFig,'name','hFig','Position',[-951,129,944,988])

%% Run plotWH
plotWH(forces,W,H,0)
figHandles = get(groot, 'Children');
set(figHandles(contains({figHandles(:).Name},{'SynFig'})),'Position',[9,9,944,1108])
set(figHandles(contains({figHandles(:).Name},{'CFig'})),'Position',[969,9,944,1108])
close(figHandles(contains({figHandles(:).Name},{'RecombFig'})))