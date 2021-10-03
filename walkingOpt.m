tstart = tic;
ds = importdata('G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\JointMotion_normalwalking.txt');
normalWalking = ds.data(2041:4108,2:4);
temp = normalWalking(:,2);
normalWalking(:,2) = normalWalking(:,3);
normalWalking(:,3) = temp;

normalWalking(:,3) = .7582.*normalWalking(:,3);

numSamples = 25;
qPos = normalWalking(floor(linspace(1,length(normalWalking),numSamples)),:);
%qPos = normalWalking([259 775 1292 1809],:);

outStim = zeros(length(qPos),7);
fVals = zeros(length(qPos),1);
outFulls = cell(length(qPos),1);
prevPos = randi(20,[1 7]);
 for ii = 1:size(qPos,1)
    tstartIn = tic;
    [outI,fVal,outFull] = legPositionOpt(qPos(ii,:),prevPos);
    outStim(ii,:) = outI;
    fVals(ii) = fVal;
    outFulls{ii} = outFull;
    prevPos = outI;
    telapsedIn = toc(tstartIn);
    disp([num2str(ii),' ',num2str(telapsedIn),' s']);
 end

telapsed = toc(tstart);
disp(['Total Calculation time:',' ',num2str(telapsed),' s'])
 
if 1
    clear outWave
    for ii = 1:length(outFulls)
        outWave(ii,:) = outFulls{ii}(end,:);
    end
    outWave = cell2mat(outWave);
    figure;
    subplot(3,1,1);plot(outWave)
    subplot(3,1,2); plot(outStim)
    subplot(3,1,3); plot(fVals)
end

% return
% outStim = [15.8129 18.4102 2.7185 14.5938 2.9031 0 1.5673;...
%             4.3689    2.1225   16.1874    7.9044    3.2207    2.6843    0.0000;...
%            4.034 9.4158 7.8383 7.5457 1.3015 6.0639 7.1284;...
%            5.3748   19.1509   17.6708    7.8534   12.8833    1.3014   10.3479;...
%            5.5444 7.0619 7.9185 7.3198 8.7896 1.6342 1.7966;...
%            4.0976   13.9300   14.7277    6.8584   10.5550    2.0589   14.6947;...
%            .2152 20 .0366 8.0055 .3919 .2388 .0064;...
%            1.1570   18.0857    7.7714    7.5301    0.4140    4.8178    6.2406;...
%            2.6257 18.3088 20 6.7828 3.6559 5.1038 .7286];
outStim3 = [outStim;outStim;outStim];
time = ds.data(:,1);

for jj = 1:size(outStim3,2)
    outStim3Big(:,jj) = interp1(1:size(outStim3,1),outStim3(:,jj),linspace(1,size(outStim3,1),length(time)));
    [fitresult] = sumsinesFit(time, outStim3Big(:,jj));
    % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
    coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
    % Equations are in the format necessary for integration into Animatlab's .asim filetype
    equations{jj} = sum_of_sines_maker(coeffs,0);
    equations_proj{jj} = sum_of_sines_maker(coeffs,1);
end
%%
simPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking_reduced_Standalone.asim";
inText = importdata(simPath);
stimInds = find(contains(inText,'<CurrentBurstOff>'));
stimNames = inText(stimInds-4);
for kk = 1:length(stimInds)
    if ~contains(inText{stimInds(kk)+1},'<CurrentOnEquation>')
        newLine = ['<CurrentOnEquation>',equations{kk},'</CurrentOnEquation>'];
        inText = [inText(1:stimInds(kk));{newLine};inText(stimInds(kk)+1:end)];
    end
end

txtInd = find(contains(inText,'OutputFilename>Joint'));
inText{txtInd} = replaceBetween(inText{txtInd},'ion','.txt','_opt1');

testSimPath = 'G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\JointMotion_opt1.asim';
fileID = fopen(testSimPath,'w');
fprintf(fileID,'%s\n',inText{:});
fclose(fileID);
sour_folder = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin';
%executable = ['"',sour_folder,'\AnimatSimulator" "',testSimPath,'"'];
executable = [string([sour_folder,'\AnimatSimulator']),string(testSimPath)];
jsystem(executable);
%subtract simPos from desPos
    dsTest = importdata('G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\JointMotion_opt1.txt');
    jointProfileTest = dsTest.data(:,3:5);
to_plot = 1;
if to_plot
    figure('Position',[1,396,1920,404])
    subplot(1,4,1)
    plot(linspace(0,3.333,length(normalWalking)),normalWalking,'LineWidth',3')
    title('Normal Walking')
    legend({'Ankle';'Knee';'Hip'})
    subplot(1,4,2)
    plot(linspace(0,3.333,length(outStim)),outStim,'LineWidth',3)
    title('Downsampled Stimulus')
    subplot(1,4,3)
    for ii = 1:7
        bb = ezplot(equations_proj{ii},[0 10/3]);
        ylim([0 20])
        bb.LineWidth = 3;
        hold on
    end
    title('Recreated stimulus waveforms')
    subplot(1,4,4)
    plot(linspace(0,3.333,length(jointProfileTest)),jointProfileTest,'LineWidth',3)
    legend({'Ankle';'Knee';'Hip'})
    title('3 Simulated Step Using Stim Waveforms')
end