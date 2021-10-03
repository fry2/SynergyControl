function h = temp_LM2020(obj,Kp)
    Ks = obj.musc_obj{2}.Kse; 
    %Kp = obj.musc_obj{2}.Kpe; 
    Lr = obj.musc_obj{2}.RestingLength; 
    STmax = obj.musc_obj{2}.ST_max;
    Fmax = obj.musc_obj{2}.max_force;
    numSamp = 100;

    [L,Am] = meshgrid(linspace(.5*Lr,1.5*Lr,numSamp),linspace(0,STmax,numSamp));
%     T = max((Ks/(Ks+Kp)).*(Kp.*(L-Lr)+(1-(L-Lr).^2./(.5*Lr)^2).*Am)./Fmax,0);
    T = (Ks/(Ks+Kp)).*(max(Kp.*(L-Lr),0)+(1-(L-Lr).^2./(.5*Lr)^2).*Am)./Fmax;
    T(T<0) = NaN;

    labelFSize = 14;

    h = figure('Position',[769.8000000000001,1.8,766.4000000000002,780.8000000000001]);
    s1 = surf(L,Am,T);
    %pbaspect([1 1 1])
    xlabel('^L/_{L_{rest}}','FontSize',labelFSize,'Position',[0.046447769721012,-1.62112229065127,-0.086187965129326])
    ylabel('^{A_{m}}/_{ST_{max}}','FontSize',labelFSize,'Position',[0.013928228823547,8.307091320371057,-0.041423164615942])
    zlabel('^{T^{*}}/_{F_{max}}','FontSize',labelFSize,'Rotation',0,'Position',[0.0088,16.30,0.41])
    plot1 = gca;
    normLticker = .5:.25:1.5;
    plot1.XTick = Lr.*normLticker;
    plot1.XTickLabel = mat2cell(normLticker,1,length(normLticker));
    normAmticker = 0:.25:1;
    plot1.YTick = STmax.*normAmticker;
    plot1.YTickLabel = mat2cell(normAmticker,1,length(normAmticker));
    s1.EdgeAlpha = .2;

    hold on

%     T2 = max((Ks/(Ks+Kp)).*(Kp.*(L-Lr)+(1-(L-Lr).^2./(.5*Lr)^2).*Am)./Fmax,0);
    T2 = (Ks/(Ks+Kp)).*(max(Kp.*(L-Lr),0)+(1-(L-Lr).^2./(.5*Lr)^2).*Am)./Fmax;      
    T2(T2>0)=NaN;
    s2 = surf(L,Am,T2);
    s2.FaceColor = 'red';
    s2.EdgeAlpha = .2;

    L3 = @(Am) Lr+(Kp*Lr^2)./(8.*Am);
    Am3 = linspace(0,STmax,numSamp);
    lLine = L3(Am3);
    lLog = lLine<=(1.6*Lr);
    lLine = lLine(lLog);
    Am3 = Am3(lLog);
    T3 = (Ks/(Ks+Kp)).*(Kp.*(lLine-Lr)+(1-(lLine-Lr).^2./(.5*Lr)^2).*Am3)./Fmax;
    plot3(lLine,Am3,T3,'c','LineWidth',5)

    plot3(Lr,STmax,1,'r.','MarkerSize',20)
    xlim([.5*Lr 1.5*Lr])
    zlim([0 1.1])
    title(num2str(Kp/obj.musc_obj{2}.Kpe))
    %view([-38 38])
    view([0 0])

    L4 = @(Am) 1+(1-sqrt(1+16*(Am./(Lr*Kp)).^2))./(8.*Am./(Lr*Kp));
    bb = L4(Am3);
end