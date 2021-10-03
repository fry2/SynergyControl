prilPath = "G:\My Drive\Rat\Swing Paper\Prilutsky15\PrilutskyCatEMGData.xlsx";
prilData = normalize_excel_data(prilPath,37);
prilData.name = 'Pril Cat';

prochPath = "G:\My Drive\Rat\Swing Paper\Prochazka\ProchazkaCatEMGData.xlsx";
prochData = normalize_excel_data(prochPath,33);
prochData.name = 'Proch Cat';

schillPath = "G:\My Drive\Rat\Swing Paper\Schilling\SchillingDogEMGData.xlsx";
schillData = normalize_excel_data(schillPath,50);
schillData.name = 'Schill Dog';

schuPath = "G:\My Drive\Rat\Swing Paper\Schubert\SchuMouseEMGData.xlsx";
schuData = normalize_excel_data(schuPath,33);
schuData.name = 'Schu Mouse';

allEMGData = struct();
allEMGData.prilData = prilData;
allEMGData.prochData = prochData;
allEMGData.schillData = schillData;
allEMGData.schuData = schuData;

%% plot data as bar
mnum = 6; expStruct = schillData; 
try
    data = expStruct.data(:,mnum); 
catch
    warning('Mnum outside range.')
return
end
[x,y] = meshgrid(1:2,linspace(0,length(data),100));
figure('Position',[1356,880,560,111]);surf(x,y,[data,data],'EdgeAlpha',0);xlim([1,2]);set(gca,'xtick',[]);view([90 90]);title([expStruct.name,' ',expStruct.cols{mnum}]);

%% Display muscle info
allDatFields = fields(allEMGData); emgDataMap = cell(1,1);
for ii = 1:length(allDatFields) % for each paper, process the muscles
    data = allEMGData.(allDatFields{ii});
    if ii == 1 % if you're the first paper, include all the muscles
        emgDataMap = data.cols';
        emgDataMap(:,2) = num2cell(1);
    else % if not, we need to check which new muscles we have
        for jj = 1:length(data.cols)
            dMusc = data.cols{jj};
            if ~any(strcmp(emgDataMap(:,1),dMusc)) % this muscle doesn't exist on the map
                emgDataMap{end+1,1} = dMusc;
                emgDataMap{size(emgDataMap,1),ii+1} = 1;
            else
                dRow = find(strcmp(emgDataMap(:,1),dMusc));
                emgDataMap{dRow,ii+1} = 1;
            end
        end
    end
end
%
[r,c] = size(emgDataMap);
for ii = 1:r
    for jj = 2:c
        if isempty(emgDataMap{ii,jj})
            emgDataMap{ii,jj} = 0;
        end
    end
    emgDataMap{ii,c+1} = sum(cell2mat(emgDataMap(ii,2:end)));
end
%
emgDataMap = sortrows(emgDataMap,c+1,'descend');
[r,c] = size(emgDataMap); temp = cell(r+1,c);
temp(2:r+1,:) = emgDataMap;
temp(1,2:5) = allDatFields';
emgDataMap = temp;
%% function import data from all sheets
function outData = process_excel_data(inSimPath)
    [~,sheets] = xlsfinfo(inSimPath);
    outData = cell(length(sheets),1);
    for ii = 1:length(sheets)
        outData{ii,1} = sheets{ii};
        outData{ii,2} = table2array(readtable(inSimPath,'Sheet',sheets{ii}));
    end
end
%% function: normalize data
function outData = normalize_excel_data(inSimPath,transInd)
    rawData = process_excel_data(inSimPath); resizedData = zeros(100,size(rawData,1));
    for ii = 1:size(rawData,1)
        temp = interp1(rawData{ii,2}(:,1),rawData{ii,2}(:,2),linspace(1,max(rawData{ii,2}(:,1)),100))';
        outRaw(:,ii) = temp';
        temp(temp < 0) = 0;
        temp = temp./max(temp);
        swing = interp1(1:transInd,temp(1:transInd),linspace(1,transInd,37));
        stance = interp1(transInd+1:length(temp),temp(transInd+1:end),linspace(transInd+1,100,63));
        resizedData(:,ii) = [swing';stance'];
    end
    outData = struct();
    outData.data = resizedData;
    outData.cols = rawData(:,1)';
    outData.rawdata = outRaw;
end