function [r2scores,recompiled,W,H] = NMFdecomposition(k,forces,to_plot,Wth)
    % Decompose an input array into component W and H arrays.
    % Input: k: integer, number of synergies to generate
    % Input: forces: array to decompose, generally forces over time
    % Input: to_plot: boolean for whether to generate plots with plotWH
    % Output: r2scores: array of r^2 values comparing the recombined waveforms with the input waveforms
    % Output: recompiled: array of same dimensions as input forces representing the array that is recombined from synergies
    % Output: W: spatial array of synergies representing muscle weighting
    % Output: H: timing array of synergies representing activation over time
    
    if size(forces,1)>size(forces,2)
        forces = forces';
    end
    
    if nargin < 4
        Wth = 0;
    end
    
    isforce = abs(max(forces,[],'all')) > 1;
    
    if isforce
        forces(forces<0)=0;
    end
    m = size(forces,1);
    %k = 4;
    gif_plot = 0;
    
    rng(300)
    W = rand(m,k);
    %H = rand(k,n);
    lamH = .01;
    lamW = .01;
    %% If conditions for plotting GIFs
    if gif_plot
        close all
        activationfig = figure('color','white','Position',[100 500 700 500]);
        weightingfig = figure('color','white','Position',[900 500 800 500]);
        recompfig = figure('color','white','Position',[50 50 1000 600]);
        count = 0;
        act_filename = [pwd,'\OutputFigures\Gifs\',datestr(datetime('now'),'yyyymmdd'),'_','cprocess.gif'];
        wgt_filename = [pwd,'\OutputFigures\Gifs\',datestr(datetime('now'),'yyyymmdd'),'_','wgtprocess.gif'];
        recom_filename = [pwd,'\OutputFigures\Gifs\',datestr(datetime('now'),'yyyymmdd'),'_','recprocess.gif'];
    end
    %% Main loop
    for i = 1:10e2
        % Multiplicative update algorithm for NMF in Berry 2006, "Algorithms and applications..."
% %         H = H.*(W'*forces)./(W'*W*H+10^-9);
% %         holder = H*H';
% %         W = W.*(forces*H')./(W*holder+10^-9);
        % ACLS algorithm for NMF from Langville 2014, "Algorithms , Initializations,..."
        % forces bit better than multiplicative update because it results in a less noisy H array
        % Additional degrees of freedom w lambda values, which affect sparsity of W and H. Ultimately little effect for this work
        % since the input arrays are not as sparse as, say, arrays for speech detection. Generally, keep lambdas [0,1]
            H = linsolve(W'*W+lamH*eye(k),W'*forces);
            H(H<0) = 0;
            W = linsolve(H*H'+lamW*eye(k),H*forces')';
            W(W./max(W)<Wth) = 0;
            % Sometimes, the solution to the local minimum equations will result in synergies being unused. To prevent this, we "bump" the offending synergy
            % with a random vector and continue the minimization
            if any(mean(W)==0)
                % find the W columns that are zero
                zW = find(mean(W)==0);
                % for each one, fill it with a random vector
                for jj = 1:length(zW)
                    W(:,zW(jj)) = rand(m,1);
                end
            end
        %% Gif addition within Loop
        if mod(i,1) == 0 && gif_plot && i < 10
            count = count + 1;
            relW = W./max(W);
            bigH = H'.*max(W);
            relW = W;
            bigH = H';
            figure(activationfig)
%             plot(linspace(0,100,size(H,2)),H','LineWidth',2)
            plot(linspace(0,100,size(H,2)),bigH,'LineWidth',2)
            if count == 1
                %ymax_act = 1.1*max(max(bigH));
                yLims = [min(bigH,[],'all') max(bigH,[],'all')];
            end
            ylim(yLims)
            drawnow
            title({'Synergy Activation Level';['NNMF Iteration: ',num2str(i)]},'FontSize',16)
            
            clf(weightingfig)
            figure(weightingfig)
            ha = tight_subplot(5,1,.05,.12,.1);
            for k = 1:size(W,2)
                axes(ha(k));
                holder = relW(:,k);
                %subplot(size(W,2),1,k)
                bar(holder)
                if k == 1
                    title({'Relative Activation of Individual Muscles';['NNMF Iteration: ',num2str(i)]},'FontSize',16)
%                     ymax = 1.1*max(W);
                elseif k == 5
                    xlabel('Muscle #')
                end
                xlim([0 39])
%                 ylim([0 ymax(k)])
                ylim([0 1.1])
                ylabel('Activation')
            end
                
            clf(recompfig)
            figure(recompfig)
            subplot(2,1,1)
            plot(forces')
            subplot(2,1,2)
            recompiled = (relW*bigH')';
            plot(recompiled)
            
            for j = 1:3
               switch j
                   case 1
                       fighand = activationfig;
                       filename = act_filename;
                   case 2
                       fighand = weightingfig;
                       filename = wgt_filename;
                   case 3
                       fighand = recompfig;
                       filename = recom_filename;
               end
              % Capture the plot as an image 
              frame = getframe(fighand); 
              im = frame2im(frame); 
              [imind,cm] = rgb2ind(im,256); 
              % Write to the GIF File 
              if count == 1 
                  imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
              else 
                  imwrite(imind,cm,filename,'gif','WriteMode','append'); 
              end
            end
        end
    end
    %% Generate r2 scores
    % r^2 value representing Pearson correlation factor
    r2scores = zeros(size(W,1),1);
    % VAF anonymous function as described by Torres-Olviedo 2006, based on uncentered Pearson correlation coefficient
    % Seems to be much more lenient about correlation verification than the r^2 values
    VAF = @(x,y) (1/length(x))*sum((x./sqrt((1/length(x))*sum(x.^2))).*(y./sqrt((1/length(y))*sum(y.^2))));
    recompiled = (W*H)';
    for i = 1:size(W,1)
        X = recompiled(:,i);
        Y = forces(i,:)';
        rho = corr(X,Y,'Type','Pearson')^2;
        %rho = VAF(X,Y);
        r2scores(i,1) = rho;
    end
    
    if any(mean(W)==0)
        keyboard
    end
    
    if to_plot
        if Wth ~= 0
            plotWH(forces,{W,Wth},H,0);
        else
            plotWH(forces,W,H,0);
        end
    end
end