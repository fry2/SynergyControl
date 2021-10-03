inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim";
%inSimPath = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_pared_Standalone.asim";

%% Test effect of length on joint VE parameters (K_elas, K_grav, and Damping)
% Some data for this has been saved in \Data since runs take a long time
    % load([pwd,'\Data\lengthVjointK.mat'],'lengthVjointK')
    % elasK = lengthVjointK.elasK; gravK = lengthVjointK.gravK; bRot = lengthVjointK.bRot; bRot_approx = lengthVjointK.bRot_approx;
    % elasK_approx = lengthVjointK.elasK_approx; gravK_approx = lengthVjointK.gravK_approx; moDoc = lengthVjointK.moDoc;
    % jmDoc = lengthVjointK.jmDoc; lengthVals = lengthVjointK.lengthVals;
lengthVals = exp(linspace(log(.6),log(100),5));
moMean = @(mo) mean(abs(diff(mo(sum(mo(1:38,:),2)~=0,100:end-100)')./.54e-3),'all');
count = 1;
elasK = zeros(18529,3,length(lengthVals)); gravK = elasK; elasK_approx = elasK; gravK_approx = elasK; bRot = elasK; bRot_approx = elasK;
jmDoc = cell(1,length(lengthVals)); moDoc = zeros(length(lengthVals),3);
for ii = 1:length(lengthVals)
    lengthScale = lengthVals(ii);
    outPath = legScaler(inSimPath,lengthScale);
    [telasK,tgravK,tbRot,zeta,telasK_approx,tgravK_approx,tbRot_approx,obj,jm] = muscle_to_joint_VE_all(outPath);
    elasK(:,:,ii) = telasK; gravK(:,:,ii) = tgravK; bRot(:,:,ii) = tbRot;
    elasK_approx(:,:,ii) = telasK_approx; gravK_approx(:,:,ii) = tgravK_approx; bRot_approx(:,:,ii) = tbRot_approx;
    disp(['Completed ',num2str(count),' out of ',num2str(length(lengthVals))])
    jmDoc{count} = jm; %moDoc(count,:) = [moMean(moH),moMean(moK),moMean(moA)];
    count = count +1;
end
baseRatLenScale = .2;
lengthVals_inM = baseRatLenScale.*lengthVals;
% Save
% lengthVjointK = struct();
% lengthVjointK.elasK = elasK; lengthVjointK.gravK = gravK; lengthVjointK.lengthVals = lengthVals; lengthVjointK.jmDoc = jmDoc;
% lengthVjointK.bRot = bRot;
% lengthVjointK.elasK_approx = elasK_approx; lengthVjointK.gravK_approx = gravK_approx; lengthVjointK.bRot_approx = bRot_approx; lengthVjointK.moDoc = moDoc;
% save([pwd,'\Data\lengthVjointK.mat'],'lengthVjointK')
%% Plot Log Scale Kvals v Length
baseRatLenScale = .2;
lengthVals_inM = baseRatLenScale.*lengthVals;
figure('Position',[962,2,958,994]);
to_plot = 0;
switch to_plot
    case 0
        ktot = elasK(100:end-100,:,:)+gravK(100:end-100,:,:);
        k2p_mean = squeeze(mean(ktot))';
         k2p_max = squeeze(max(ktot,[],1))';
         k2p_min = squeeze(min(ktot,[],1))';
%        k2p_max = k2p_mean+squeeze(std(ktot));
%        k2p_min = k2p_mean-squeeze(std(ktot));
        titlePhrase = 'Joint Total Stiffness v. Length Scale';
        legPhrase = {'Mean';'Max';'Min'};
    case 1
        k2p_mean = squeeze(mean(elasK(100:end-100,:,:)));
        k2p_max = squeeze(max(elasK(100:end-100,:,:),[],1));
        k2p_min = squeeze(min(elasK(100:end-100,:,:),[],1));
        titlePhrase = 'Joint Elastic Stiffness v. Length Scale';
        legPhrase = {'Mean';'Max';'Min'};
    case 2
        k2p_mean = squeeze(mean(gravK(100:end-100,:,:)));
        k2p_max = squeeze(max(gravK(100:end-100,:,:),[],1));
        k2p_min = squeeze(min(gravK(100:end-100,:,:),[],1));
        titlePhrase = 'Joint Gravitational Stiffness v. Length Scale';
        legPhrase = {'Mean';'Max';'Min'};
    case 3
        ktot = elasK(100:end-100,:,:)+gravK(100:end-100,:,:);
        krat2 = elasK(100:end-100,:,:)./ktot;
        krat1 = gravK(100:end-100,:,:)./ktot;
        k2p_mean1 = squeeze(mean(krat1));
        k2p_mean2 = squeeze(mean(krat2));
%         k2p_max = squeeze(max(krat1,[],1));
%         k2p_min = squeeze(min(krat1,[],1));
        titlePhrase = 'Proportion of Elastic Stiffness Over Total Stiffness v. Length Scale';
        legPhrase = {'k_{grav}/k_{tot}';'k_{elas}/k_{tot}'};
    case 4
        k2p_mean = squeeze(mean(bRot(100:end-100,:,:)));
        k2p_max = squeeze(max(bRot(100:end-100,:,:),[],1));
        k2p_min = squeeze(min(bRot(100:end-100,:,:),[],1));
        titlePhrase = 'Joint Damping v. Length Scale';
        legPhrase = {'Mean';'Max';'Min'};
end
subplot(3,1,1)
    if to_plot == 3
        semilogx(lengthVals_inM,k2p_mean1(:,1),'LineWidth',3);hold on;grid on;
        semilogx(lengthVals_inM,k2p_mean2(:,1),'LineWidth',3);
    else
        loglog(lengthVals_inM,k2p_mean(:,1),'LineWidth',3);hold on;grid on;
        loglog(lengthVals_inM,k2p_max(:,1));
        loglog(lengthVals_inM,k2p_min(:,1));
    end
    title(titlePhrase,'FontSize',16)
    ylabel('Hip Stiffness (N-m)','FontSize',16)
    xlabel('Length Scale (m)','FontSize',16)
    xlim([min(lengthVals_inM) max(lengthVals_inM)])
    legend(legPhrase,'Location','northwest')
subplot(3,1,2)
    if to_plot == 3
        semilogx(lengthVals_inM,k2p_mean1(:,2),'LineWidth',3);hold on;grid on;
        semilogx(lengthVals_inM,k2p_mean2(:,2),'LineWidth',3);
    else
        loglog(lengthVals_inM,k2p_mean(:,2),'LineWidth',3);hold on;grid on;
        loglog(lengthVals_inM,k2p_max(:,2));
        loglog(lengthVals_inM,k2p_min(:,2));
    end
    ylabel('Knee Stiffness (N-m)','FontSize',16)
    xlabel('Length Scale (m)','FontSize',16)
    xlim([min(lengthVals_inM) max(lengthVals_inM)])
    legend(legPhrase,'Location','northwest')
subplot(3,1,3)
    if to_plot == 3
        semilogx(lengthVals_inM,k2p_mean1(:,3),'LineWidth',3);hold on;grid on;
        semilogx(lengthVals_inM,k2p_mean2(:,3),'LineWidth',3);
    else
        loglog(lengthVals_inM,k2p_mean(:,3),'LineWidth',3);hold on;grid on;
        loglog(lengthVals_inM,k2p_max(:,3));
        loglog(lengthVals_inM,k2p_min(:,3));
    end
    ylabel('Ankle Stiffness (N-m)','FontSize',16)
    xlabel('Length Scale (m)','FontSize',16)
    xlim([min(lengthVals_inM) max(lengthVals_inM)])
    legend(legPhrase,'Location','northwest')
%% Plot Joint VE parameters v Joint Angle
figure('Position',[962,2,958,994]);
only_rat = 1;
if only_rat
    ratInd = find(lengthVals==1+min(abs(lengthVals-1)));
end
to_plot = 0;
switch to_plot
    case 0
       ve2plot = elasK(100:end-100,:,:)+gravK(100:end-100,:,:);
       ve2plotA = elasK_approx(100:end-100,:,:)+gravK_approx(100:end-100,:,:);
       titlePhrase = 'Total Joint Stiffness v Joint Angle';
       ylabelPhrase = ' Joint Stiffness (N-m)';
    case 1
       ve2plot = elasK(100:end-100,:,:);
       ve2plotA = elasK_approx(100:end-100,:,:);
       titlePhrase = 'Elastic Joint Stiffness v Joint Angle';
       ylabelPhrase = ' Joint Stiffness (N-m)';
    case 2
       ve2plot = gravK(100:end-100,:,:);
       ve2plotA = gravK_approx(100:end-100,:,:);
       titlePhrase = 'Gravitational Joint Stiffness v Joint Angle';
       ylabelPhrase = ' Joint Stiffness (N-m)';
    case 3 
       ve2plot = bRot(100:end-100,:,:);
       ve2plotA = bRot_approx(100:end-100,:,:);
       titlePhrase = 'Joint Damping v Joint Angle';
       ylabelPhrase = ' Joint Damping (Ns-m)';
end
subplot(3,1,1)
[a,b] = size(ve2plot);
    switch only_rat
        case 1
            plot(jmDoc{ratInd}(100:end-100,1),1000.*ve2plotA(:,ratInd,1),'r','LineWidth',1); grid on; hold on
            plot(jmDoc{ratInd}(100:end-100,1),1000.*ve2plot(:,ratInd,1),'b','LineWidth',4); 
            legend({'Small Angle';'Exact'},'Location','northwest')
            if to_plot ~= 3
                ylabel(['Hip',ylabelPhrase(1:end-5),'(mNm)'],'FontSize',16)
            else
                ylabel(['Hip',ylabelPhrase],'FontSize',16)
            end
        otherwise
            for ii = 1:b
                plot(jmDoc{ii}(100:end-100,1),ve2plot(:,ii,1)); hold on
            end
            ylabel(['Hip',ylabelPhrase],'FontSize',16)
    end
    title(titlePhrase,'FontSize',16)
    xlabel('Joint Angle (deg)')
    xlim([min(jmDoc{5}(:,1)) max(jmDoc{5}(:,1))])
subplot(3,1,2)
    switch only_rat
        case 1
            plot(jmDoc{ratInd}(100:end-100,2),1000.*ve2plotA(:,ratInd,2),'r','LineWidth',1);grid on; hold on
            plot(jmDoc{ratInd}(100:end-100,2),1000.*ve2plot(:,ratInd,2),'b','LineWidth',4); 
            legend({'Small Angle';'Exact'},'Location','northwest')
            if to_plot ~=3
                ylabel(['Knee',ylabelPhrase(1:end-5),'(mNm)'],'FontSize',16)
            else
                ylabel(['Knee',ylabelPhrase],'FontSize',16)
            end
        otherwise
            for ii = 1:b
                plot(jmDoc{ii}(100:end-100,2),ve2plot(:,ii,2)); hold on
            end
            ylabel(['Knee',ylabelPhrase],'FontSize',16)
    end
    xlabel('Joint Angle (deg)')
    xlim([min(jmDoc{5}(:,2)) max(jmDoc{5}(:,2))])
subplot(3,1,3)
    switch only_rat
        case 1
            plot(jmDoc{ratInd}(100:end-100,3),1000.*ve2plotA(:,ratInd,3),'r','LineWidth',1); grid on; hold on
            plot(jmDoc{ratInd}(100:end-100,3),1000.*ve2plot(:,ratInd,3),'b','LineWidth',4); 
            legend({'Small Angle';'Exact'},'Location','northwest')
            if to_plot ~= 3    
                ylabel(['Ankle',ylabelPhrase(1:end-5),'(mNm)'],'FontSize',16)
            else
                ylabel(['Ankle',ylabelPhrase],'FontSize',16)
            end
        otherwise
            for ii = 1:b
                plot(jmDoc{ii}(100:end-100,3),ve2plot(:,ii,3)); hold on
            end
            ylabel(['Ankle',ylabelPhrase],'FontSize',16)
    end
    xlabel('Joint Angle (deg)')
    xlim([min(jmDoc{5}(:,3)) max(jmDoc{5}(:,3))])
    clear a b ve2plot* titlePhrase ratInd legPhrase to_plot ylabelPhrase
%% Plot S Vals v Joint Angles
% S is the ratio of moment arm over muscle length. Sutton sets this to a constant .032, rather than a dynamic value
musc2check = sum(moH(1:38,:),2)~=0;
hipS = abs(moH(musc2check,:)')./mlH(:,musc2check);
musc2check = sum(moK(1:38,:),2)~=0;
kneeS = abs(moK(musc2check,:)')./mlK(:,musc2check);
musc2check = sum(moA(1:38,:),2)~=0;
ankleS = abs(moA(musc2check,:)')./mlA(:,musc2check);

% S relative to bone length
% musc2check = sum(moH(1:38,:),2)~=0;
% hipS = abs(moH(musc2check,:)')./bodyLengths(1);
% musc2check = sum(moK(1:38,:),2)~=0;
% kneeS = abs(moK(musc2check,:)')./bodyLengths(1);
% musc2check = sum(moA(1:38,:),2)~=0;
% ankleS = abs(moA(musc2check,:)')./bodyLengths(2);
    
figure('Position',[962,2,958,994]);
    subplot(3,1,1)
        plot(jm(100:end-100,1),hipS(100:end-100,:),'LineWidth',1.5); grid on; hold on
        plot(jm(100:end-100,1),.032.*ones(length(jm(100:end-100,1)),1),'r','LineWidth',3)
        title('Moment Arm/Muscle Length v. Joint Angle','FontSize',18)
        ylabel('r/L','FontSize',16)
        xlabel('Hip Joint Angle (deg)','FontSize',16)
        xlim([min(jm(100:end-100,1)), max(jm(100:end-100,1))])
    subplot(3,1,2)
        plot(jm(100:end-100,2),kneeS(100:end-100,:),'LineWidth',1.5); grid on; hold on
        plot(jm(100:end-100,2),.032.*ones(length(jm(100:end-100,1)),1),'r','LineWidth',3)
        ylabel('r/L','FontSize',16)
        xlabel('Knee Joint Angle (deg)','FontSize',16)
        xlim([min(jm(100:end-100,2)), max(jm(100:end-100,2))])
    subplot(3,1,3)
        plot(jm(100:end-100,3),ankleS(100:end-100,:),'LineWidth',1.5); grid on; hold on
        plot(jm(100:end-100,3),.032.*ones(length(jm(100:end-100,1)),1),'r','LineWidth',3)
        ylabel('r/L','FontSize',16)
        xlabel('Ankle Joint Angle (deg)','FontSize',16)
        xlim([min(jm(100:end-100,3)), max(jm(100:end-100,3))])
%% Plot Zeta for Each Joint
    bodyMasses = [14.141,3.342,1.571].*10^-3;
    bodyLengths = [35.7,44.9,21.3].*10^-3; % body lengths in meters (Fe,T,Fo)
    m = zeros(1,3);
    for ii = 1:3
        m(ii) = (1/3)*bodyLengths(ii)^2*sum(bodyMasses(ii:end));
    end
    [temph,hjm,tempha,moH,mlH] = muscle_to_joint_VE(inSimPath,1);
    [tempk,kjm,tempka,moK,mlK] = muscle_to_joint_VE(inSimPath,2);
    [tempa,ajm,tempaa,moA,mlA] = muscle_to_joint_VE(inSimPath,3);
    
    c = [temph(:,2),tempk(:,2),tempa(:,2)]; k = [temph(:,1)+temph(:,3),tempk(:,1)+tempk(:,3),tempa(:,1)+tempa(:,3)];
    
    zeta = c./(2.*sqrt(m.*k));

    for ii = 1:40
        m = zeros(1,3);
        for jj = 1:3
            m(jj) = (1/3)*bodyLengths(jj)^2*sum(bodyMasses(jj:end)).*lengthVals(ii).^3;
        end
        c = squeeze(bRot(:,ii,:)); k = squeeze(elasK(:,ii,:)+gravK(:,ii,:));
        zeta = c./(2.*sqrt(m.*k));
        zetaScale(ii,:) = mean(zeta);
    end
    
    figure;semilogx(lengthVals_inM,zetaScale,'LineWidth',3); xlabel('Length (m)'); ylabel('Zeta')