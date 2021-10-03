close all
Kp1 = obj.musc_obj{2}.Kpe;

%temp_LM2020(obj,Kp);
%[Tension,Tdot,Act,mL] = animatlab_CalculateTension(obj,obj.musc_obj{2},Kp);

for ii = 0:.1:1.5
    Kp = ii*Kp1;
    h = temp_LM2020(obj,Kp);
    %animatlab_CalculateTension(obj,obj.musc_obj{2},Kp);
    %h = gcf;
    %saveas(h,['G:\My Drive\Rat\SynergyControl\OutputFigures\Images\Kp_change2\',erase(num2str(ii),'.'),'.png']);
end