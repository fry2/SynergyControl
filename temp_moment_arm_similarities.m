% assumes you already have a FullLeg obj
 sim_file = 'C:\Users\fry16\OneDrive\Documents\JointDampingOpt\JointDampingOpt_compAnkle_Standalone.asim';
 obj = FullLeg(sim_file,[],[]);
%%
nicholslist = {'lateralgastrocnemius';...
'plantaris';...
'medialgastrocnemius';...
'flexorhallucislongus';...
'soleus';...
'flexordigitorumlongus';...
'tibialisposterior';...
'pectineus';...
'bicepsfemorisanterior';...
'semimembranosus';...
'adductorlongus';...
'bicepsfemorisposterior';...
'gracilisanterior';...
'semitendinosusprincipal';...
'quadratusfemoris';...
'gluteusmaximus';...
'gluteusmedius';...
'gluteusminimus';...
'extensordigitorumlongus';...
'tibialisanterior';...
'peroneuslongus';...
'rectusfemoris';...
'vastusintermedius';...
'vastuslateralis';...
'vastusmedialis';...
'illiopsoas';...
'gracilisposterior';...
'tensorfascialatae';...
'obturatorexternus';...
'piriformis';...
'gemellussuperior';...
'caudofemoralis';...
'popliteus';...
'adductormagnus';...
'adductorbrevis';...
'gemellusinferior';...
'obturatorinternus';
'semitendinosusaccessory'};

%% Adapted corrFig
    % Animatlab Order
    mOut = zeros(38,3);
    mOutNormal = zeros(38,3);
    colcount = 1;
    for joint = 1:3
        for axis1 = 1:3
            moment_output = compute_joint_moment_arms(obj,joint,axis1);
            mOutNormal(:,colcount) = moment_output(1:38,1);
            colcount = colcount + 1;
        end
    end
    mOut = mOutNormal;
%%
    order = 1;
    mNamesNormal = cell(38,1);
    % Reorder moment arms based on preferences
    switch order
        case 0
            % Normal order, don't change anything
            matOrder = 1:38;
            for ii = 1:38
                mNamesNormal{ii} = obj.musc_obj{ii}.muscle_name(4:end);
            end
        case 1
            % Order the same as Nichols paper figure (for as many matching muscles as is possible)
            for ii = 1:38
                mNamesNormal{ii} = obj.musc_obj{ii}.muscle_name(4:end);
                nicholsOrder(ii) = find(contains(nicholslist,mNamesNormal{ii}));
            end
            matOrder = nicholsOrder;
            %mNames = nicholslist;
        case 2
            Z = linkage(mOutNormal,'weighted');[H,T,outperm] = dendrogram(Z,38,'Orientation','left','ColorThreshold',0.25*max(Z(:,3)));
            set(H,'LineWidth',2)
            set(gcf,'Position',[13,86,404,828])
            set(gca, 'YDir','reverse')
            for ii = 1:38
                mNamesNormal{ii} = obj.musc_obj{ii}.muscle_name(4:end);
                matOrder(outperm(ii)) = ii;
            end
    end

    temp = {};
    for ii = 1:38 
        mOut(matOrder(ii),:) = mOutNormal(ii,:);
        temp{matOrder(ii),1} = mNamesNormal{ii};
    end
    mNames = temp;
    
    mMatOut = zeros(38,38);
    for musc1 = 1:38
        v1 = mOut(musc1,:);
        for musc2 = 1:38
            v2 = mOut(musc2,:);
            mMatOut(musc1,musc2) = dot(v1,v2)/(norm(v1)*norm(v2));
        end
    end
    mMatOut = normalize(mMatOut,'range',[-1 1]);
    if order == 1
        mMatOut = mMatOut(1:26,1:26);
    end
    graphedInfo = [mNames,num2cell(mOut)];
% Create figure
    corrFig = figure('name','CorrFig','Position',[392,2,958,994]);
    sorted = 0;
    
    % Create axes
    axes1 = axes('Parent',corrFig);
    hold(axes1,'on');
    clims = [min(mMatOut,[],'all') max(mMatOut,[],'all')];
    axisticker = 1:length(mMatOut);
    ylabels = num2cell(axisticker);
    ylabels = mNames;
    xlabels = ylabels;
    
    if sorted == 1
        [mMatOut,newRowInds] = sortrows(mMatOut,1,'descend');
        ylabels = num2cell(newRowInds);
    elseif sorted == 2
        [mMatOut,newRowInds] = sortrows(mMatOut,1,'descend');
        [mMatOut,newColInds] = sortrows(mMatOut',1,'descend');
        mMatOut = mMatOut';
        ylabels = num2cell(newRowInds);
        for ii = 1:length(newRowInds)
            ylabels{ii} = mNames{newRowInds(ii)};
        end
        xlabels = num2cell(newColInds);
        for ii = 1:length(newColInds)
            xlabels{ii} = mNames{newColInds(ii)};
        end
    end
    
    % Create image
    image(mMatOut,'Parent',axes1,'CDataMapping','scaled');
    title('Data Correlation','FontSize',15)
   
    % Uncomment the following line to preserve the X-limits of the axes
    % xlim(axes1,[0.5 10.5]);
    % Uncomment the following line to preserve the Y-limits of the axes
    % ylim(axes1,[0.5 10.5]);
    box(axes1,'on');
    axis(axes1,'ij');
    % Set the remaining axes properties
    xlim([.5 size(mMatOut,2)+.5])
    ylim([.5 size(mMatOut,1)+.5])
    %set(gca,'xaxisLocation','top')
    pbaspect([38 38 1])
    set(axes1,'CLim',clims,'Colormap',turbo(40),'Layer','top','XTick',axisticker,'XTickLabelRotation',90,'YTick',axisticker,'XTickLabels',xlabels,'YTickLabels',ylabels);
    colorbar