clear meandiffs
count = 1;
wthreshs = linspace(0,.9,10);
stddiffs = zeros(length(wthreshs),1);
gif_plot = 0;

if gif_plot
    close all
    w_filename = [pwd,'\OutputFigures\Gifs\Wthresh\',datestr(datetime('now'),'yyyymmdd'),'_','wthresh.gif'];
    h_filename = [pwd,'\OutputFigures\Gifs\Wthresh\',datestr(datetime('now'),'yyyymmdd'),'_','hthresh.gif'];
    re_filename= [pwd,'\OutputFigures\Gifs\Wthresh\',datestr(datetime('now'),'yyyymmdd'),'_','recthresh.gif'];
    count2 = 0;
end

rng(300)

meandiffs = zeros(length(wthreshs),10);
synfocus = 5;
numSyns = 10;
caseNum = 1;

figHandles = get(groot, 'Children');
if ~isempty(figHandles)
    priorFigures = contains({figHandles(:).Name},{['surffig',num2str(caseNum)],'stdfig'});
    close(figHandles(priorFigures))
end
scrSz = get(groot, 'ScreenSize');
scW = scrSz(3);

switch caseNum
    case 1 % Normal injected current analysis
        inArray = current2inject';
    case 2 % Rearrange the muscle waveforms before analysis
        muscs_rand = muscs(randperm(length(muscs)));
        inArray = current2inject(muscs_rand,:)';
    case 3 % Input random noise
        inArray = rand(500,38);
end

for tt = 1:numSyns
    count = 1;
    for ii = wthreshs
        [~,recompiled,W,H] = NMFdecomposition(tt,inArray,0,ii);
        differ = abs((inArray-recompiled));
        meandiffs(count,tt) = mean(mean(differ));
        stddiffs(count,tt) = std(mean(differ));
        count = count + 1;
        if tt == synfocus && gif_plot
            count2 = count2 +1;
            if count2 > 1
                rscores = corr(Hhold',H');
                arrLog = rscores==max(max(rscores,0),[],2);
                dupcols = find(sum(arrLog)>1);
                transMat = zeros(size(arrLog));
                transMat(find(~sum(arrLog(:,dupcols),2),5),:) = arrLog(find(~sum(arrLog(:,dupcols),2),5),:);
                for jj = 1:length(dupcols)
                    maxincol = max(rscores(:,dupcols(jj)));
                    rowhold = find(rscores(:,dupcols(jj))==maxincol);
                    transMat(:,dupcols(jj)) = 0;
                    transMat(rowhold,dupcols(jj)) = 1;
                end
                emptyRows = find(~sum(transMat,2),5);
                for jj = 1:length(emptyRows)
                    sortedrs = sort(rscores(emptyRows(jj),:),'descend');
                    valtry = 1;
                    colhold = find(rscores(emptyRows(jj),:)==sortedrs(valtry));
                    while sum(transMat(:,colhold))~=0
                        colhold = find(rscores(emptyRows(jj),:)==sortedrs(valtry));
                        valtry = valtry+1;
                    end
                    transMat(emptyRows(jj),colhold) = 1;
                end
                [oldSyns,newSyns] = find(transMat,5);
                fixedH = zeros(size(H));
                fixedW = zeros(size(W));
                for dd = 1:length(transMat)
                    fixedH(oldSyns(dd),:) = H(newSyns(dd),:);
                    fixedW(:,oldSyns(dd)) = W(:,newSyns(dd));
                end
                H = fixedH;
                W = fixedW;
            end
            Hhold = H;
            plotWH(inArray,{W,ii},H,0)
            figHandles = get(groot, 'Children');
            if ~isempty(figHandles)
                w_fig = contains({figHandles(:).Name},{'SynFig'});
                h_fig = contains({figHandles(:).Name},{'CFig'});
                re_fig = contains({figHandles(:).Name},{'RecombFig'});
            end
            % Capture the plot as an image
            figs = {w_fig;h_fig;re_fig};
            figfiles = {w_filename;h_filename;re_filename};
            for jj = 1:length(figs)
                fighand = figs{jj};
                frame = getframe(figHandles(fighand)); 
                im = frame2im(frame); 
                [imind,cm] = rgb2ind(im,256); 
                % Write to the GIF File 
                if count2 == 1 
                  imwrite(imind,cm,figfiles{jj},'gif', 'DelayTime',0.8, 'Loopcount',inf); 
                else 
                  imwrite(imind,cm,figfiles{jj},'gif','DelayTime',0.8,'WriteMode','append'); 
                end
            end
        end
    end
end

    figure('name',['surffig',num2str(caseNum)],'Position',[-(scW-10)/2 130 scW/2 985])
    surf(1:numSyns,wthreshs,meandiffs)
    pbaspect([1 1 1])
    view([43 27])
    title('Differences in NNMF signals vs. Input Signals')
    xlabel('Number of Synergies')
    ylabel('Synergy Threshold Value')
    zlabel('Average Difference Between Signals')
    [minthresh,minsyn] = minmat(stddiffs);
    disp(['Minimum NNMF Error with ',num2str(minsyn),' synergies and a threshold value of ',num2str(wthreshs(minthresh)),'.'])
    
    function [ a,b ] = minmat(c)
        as=size(c);
        [~,I]=min(c(:));
        r=rem(I,as(1));
        a=r;
        b=((I-a)/as(1))+1;
        if a==0
            a=as(1);
            b=b-1;
        else
            a=r;
        end
    end
