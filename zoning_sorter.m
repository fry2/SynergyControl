function muscZones = zoning_sorter(simText,numZones)
    if ischar(simText) || isstring(simText)
        simText = importdata(simText);
    end
    mInds = find(contains(simText,'<Type>LinearHillMuscle</Type>'));
    muscZones = cell(38,2);
    for ii = 1:length(mInds)
        attdir = {};
        attsInd = find(contains(simText(mInds(ii):end),'<Attachments>'),1,'first')+mInds(ii)-1;
        attcounter = 1;
        while ~contains(simText(attsInd+attcounter),'</Attachments>') && attcounter <= 10
            attdir{attcounter,1} = char(extractBetween(string(simText{attsInd+attcounter}),'>','</'));
            attcounter = attcounter + 1;
        end
        for jj = 1:size(attdir,1)
            attName = simText{find(contains(simText,['<ID>',attdir{jj},'</ID>']))-1};
            attdir{jj,2} = find([contains(attName,'Pelvis ') contains(attName,'Femur ') contains(attName,'Tibia ') contains(attName,'Foot ')]);
            attdir{jj,3} = attName(strfind(attName,'LH_'):end-7);
        end
        muscZones{ii,1} = char(extractBetween(string(simText{mInds(ii)-2}),'>','</'));
        if numZones == 6
            if all(cellfun(@(x) isequal(x,attdir{1,3}),attdir(:,3)))
                switch attdir{1,3}
                    case 'LH_HipZ Ext'
                        muscZones{ii,2} = 1;
                    case 'LH_HipZ Flx'
                        muscZones{ii,2} = 2;
                    case 'LH_Knee Ext'
                        muscZones{ii,2} = 3;
                    case 'LH_Knee Flx'
                        muscZones{ii,2} = 4;
                    case 'LH_AnkleZ Ext'
                        muscZones{ii,2} = 5;
                    case 'LH_AnkleZ Flx'
                        muscZones{ii,2} = 6;
                end
            else
                % Some specific zoning designations for multi-purpose muscles
                if isequal(attdir(:,3),{'LH_HipZ Flx';'LH_Knee Ext';'LH_AnkleZ Ext';'LH_Knee Ext'})
                    % Rectus femoris
                    muscZones{ii,2} = 3;
                elseif isequal(attdir(:,3),{'LH_Knee Flx';'LH_Knee Ext';'LH_Knee Ext'})
                    % BFP
                    muscZones{ii,2} = 4;
                else
                    keyboard
                end
            end
        elseif numZones == 38
            muscZones{ii,2} = ii;
        else
            muscZones{ii,2} = [];
        end
    end
end