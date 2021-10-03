% Nick
c0 = 1.3e3; k0 = 12e3; m0 = 12; L0 = 1; s = sqrt(1e-3);

% Fletcher
%c0 = 57.28; k0 = 23.42; m0 = 16.94e-4; L0 = .2; s = .2;
sTemp = @(L) 10.^(.434.*log(12.*L.^3))./(100*(.1375.*L-7.5e-3));
phi = @(c_0,k_0,m_0,L,L_0,s,T) atand(((c_0.*s.^2.*(L.^3./L_0).*(2.*pi./T))./(k_0.*s.^2.*(L.^3./L_0)+(1/2).*m_0.*9.81.*(L.^4./L_0.^3)-(1/3).*m_0.*(L.^5./L_0.^3).*(2.*pi./T).^2)));
tqv = @(c_0,k_0,m_0,L,L_0,s) (pi.*c_0.*s.^2.*L.^3+pi.*((c_0.*s.^2.*L.^3).^2+16.*(k_0.*s.^2.*L.^3+(1/2).*m_0.*9.81.*(L.^4./L_0.^2)).*((1/3).*m_0.*(L.^5./L_0.^2))).^(1/2))./(2.*(k_0.*s.^2.*L.^3+(1/2).*m_0.*9.81.*(L.^4./L_0.^2)));
tvk = @(c_0,k_0,m_0,L,L_0,s) (-pi.*c_0.*s.^2.*L.^3+pi.*sqrt((c_0.*s.^2.*L.^3).^2+16.*(k_0.*s.^2.*L.^3+(1/2).*m_0.*9.81.*L.^4).*((1/3).*m_0.*(L.^5./L_0.^2))))./(2.*(k_0.*s.^2.*L.^3+(1/2).*m_0.*9.81.*(L.^4./L_0.^2)));
zta = @(c_0,k_0,m_0,L,L_0,s) (1./L).*((c_0.*s.^2)./(2.*sqrt(k_0.*s.^2+(1/2).*m_0.*9.81.*(L./L_0.^2)).*sqrt((1/3).*(m_0./L_0.^2))));

[lmat,tmat] = meshgrid(logspace(-4,1,500),logspace(-3,2,500)); 
sVals = logspace(-3,0,10);
temp = find(sVals-sqrt(1e-3)>0,1,'first');
sVals = [sVals(1:temp-1),sqrt(1e-3),sVals(temp:end)];
%sVals(1) = s;
%sVals = repmat(sTemp(logspace(-4,1,500)),500,1);

count  = 1;
for ii = 1:length(sVals)
    figObj = figure('Position',[962,2,958,994]);
    phiout = phi(c0,k0,m0,lmat,L0,sVals(ii),tmat);
    phiout(phiout<0) = phiout(phiout<0)+180;

     %surf(tmat,lmat,phiout,'EdgeAlpha',0); hold on
    pPic = pcolor(tmat,lmat,phiout); 
        pPic.EdgeAlpha = 0; hold on
        set(gca,'xscale','log','yscale','log'); 
        view([0 90]);
        pbaspect([1 1 1]);
        colormap autumn

    title(['Phase Angle of Steady State Motion, $\phi$. Sval = ',num2str(sVals(ii))],'Interpreter','latex','FontSize',18)
    xlabel('Time-scale (s)','FontSize',16);
    ylabel('Length-scale (m)','FontSize',16);

    % Plot Tqv
        plot(tqv(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii)),logspace(-4,1,50),'w','LineWidth',3)

    % Plot Tvk
        plot(tvk(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii)),logspace(-4,1,50),'w','LineWidth',3)
    
    % Plot Zeta
        zta_temp = @(L) zta(c0,k0,m0,L,L0,sVals(ii))-1;
        ztaVal = fsolve(zta_temp,.1,optimoptions('fsolve','Display','none'));
        yline(ztaVal,'k:','LineWidth',2)
        
    % Plot Rat
        plot(linspace(.2,.6,50),.05.*ones(1,50),'k','LineWidth',3)
        text(.6,.05,'\leftarrow Rat','FontSize',15)
        
    % Plot Horse
        plot(linspace(.45,.9,50),1.5.*ones(1,50),'k','LineWidth',3)
        text(.9,1.5,'\leftarrow Horse','FontSize',15)
        
    % Plot Fly
        plot(linspace(.05,.15,50),4e-4.*ones(1,50),'k','LineWidth',3)
        text(.15,4e-4,'\leftarrow Fly','FontSize',15)
    
    % Plot Aplysia
        plot(linspace(1.2,5,50),.015.*ones(1,50),'k','LineWidth',3)
        text(5,.015,'\leftarrow Aplysia','FontSize',15)
% 
    saveas(figObj,['G:\My Drive\Rat\S Analysis\Suttinski S\',num2str(count),'.png'])
    close(figObj)
    count  = count +1;
end

%% Plot Kinetic and Viscoelastic Region Sizes v. S Vals

sVals = logspace(-3,0,1000);
fullSpace = (10^2-10^-3)*(10^1-10^-4);
for ii = 1:length(sVals)
    kinSizeInit = trapz(logspace(-4,1,50),tvk(c0,k0,m0,logspace(-4,1,50),L0,sqrt(1e-3)));
    kinSize(ii) = (trapz(logspace(-4,1,50),tvk(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii)))./kinSizeInit).*100;
    vSizeInit = trapz(logspace(-4,1,50),tqv(c0,k0,m0,logspace(-4,1,50),L0,sqrt(1e-3))-tvk(c0,k0,m0,logspace(-4,1,50),L0,sqrt(1e-3)));
    vSize(ii) = ((trapz(logspace(-4,1,50),tqv(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii))-tvk(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii))))./vSizeInit).*100;
    qSizeInit = fullSpace - vSizeInit;
    qSize(ii) = ((fullSpace-(trapz(logspace(-4,1,50),tqv(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii))-tvk(c0,k0,m0,logspace(-4,1,50),L0,sVals(ii)))))./qSizeInit).*100;
end

figure('Position',[962,2,958,994]);
subplot(3,1,1)
    semilogx(sVals,kinSize,'LineWidth',3); grid on; hold on
    semilogx(.032,100,'r.','MarkerSize',30)
    title('Kinematic Region v. S Val','FontSize',18)
    ylabel('Percent of Original Size (%)')
    xlabel('S Value')
subplot(3,1,2)
    semilogx(sVals,vSize,'LineWidth',3); grid on; hold on
    semilogx(.032,100,'r.','MarkerSize',30)
    title('Viscoelastic Region v. S Val','FontSize',18)
    ylabel('Percent of Original Size (%)')
    xlabel('S Value')
subplot(3,1,3)
    semilogx(sVals,qSize,'LineWidth',3); grid on; hold on
    semilogx(.032,100,'r.','MarkerSize',30)
    title('Quasi-Static Region v. S Val','FontSize',18)
    ylabel('Percent of Original Size (%)')
    xlabel('S Value')
%% Plot damping ratio over length and Sval ranges

[sMat,lMat] = meshgrid(logspace(-3,0,500),logspace(-2,0,500));
zzMat = zta(c0,k0,m0,lMat,L0,sMat);
surf(sMat,lMat,zzMat,'EdgeAlpha',0); grid on; 
xlabel('S');    
ylabel('L');
zlabel('Zeta');
set(gca,'xscale','log','yscale','log','zscale','log');