variables = whos;
varNames = {variables(:).name};
scrSz = get(groot, 'ScreenSize');
scW = scrSz(3);

if ~any(contains(varNames,'sdata'))
    sdata = processSimData("G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim",1);
end

muscnum = 8;
dstart = 1;
mndataInd = find(ismember(sdata(3).data_cols,['neur',num2str(muscnum)]));
mnsim_results = sdata(3).data(dstart:end-10,mndataInd);
muscinfoInd = find(contains(sdata(4).data_cols,obj.musc_obj{muscnum}.muscle_name(4:end)));
muscinfo = sdata(4).data(dstart:end-10,muscinfoInd);
actinfo = sdata(6).data(dstart:end-10,muscinfoInd);
alinfo = sdata(7).data(dstart:end-10,muscinfoInd);
input_res = V_musc(muscnum,:)';

yMax = max([max(mnsim_results) max(muscinfo(10:end-10))]);
yMin = min([min(mnsim_results) min(muscinfo(10:end-10))]);

figHandles = get(groot, 'Children');
if ~isempty(figHandles)
    priorFigures = contains({figHandles(:).Name},{obj.musc_obj{muscnum}.muscle_name(4:end),'inwaves'});
    close(figHandles(priorFigures))
end

figure('Position',[-(scW-10) 130 scW/2 985],'name',obj.musc_obj{muscnum}.muscle_name(4:end))
subplot(3,1,1)
    calcv = V_musc(muscnum,1:end-20)';
    time = sdata(2).time(dstart:end-10);
    calcv = interp1(1:length(calcv),calcv,linspace(1,length(calcv),length(muscinfo)))';
    plot(time,1000*calcv,'LineWidth',2)
    hold on
    plot(time,1000*mnsim_results,'LineWidth',2)
    ylim(1000*[min([calcv;mnsim_results]) max([calcv;mnsim_results])])
    legend({'Calculated';'Simulated'},'Location','north')
    title(['MN Simulation results; ',['neur',num2str(muscnum)]])
subplot(3,1,2)
    muscle = obj.musc_obj{muscnum};
%     STmax = muscle.ST_max;
%     xoff = muscle.x_off;
%     yoff = muscle.y_off;
%     steep = muscle.steepness;
%     st_curve = @(Fmax,steepness,xoff,V,yoff) Fmax./(1+exp(steepness*(xoff-V)))+yoff;
%     sim_tension = st_curve(STmax,steep,xoff,muscinfo(10:end-10),yoff);
    sim_tension = actinfo;
    calcam = Am_musc(muscnum,1:end-10)';
    calcam = interp1(1:length(calcam),calcam,linspace(1,length(calcam),length(sim_tension)))';
    plot(calcam,'LineWidth',2)
    hold on
    plot(sim_tension,'LineWidth',2)
    legend({'Calculated';'Simulated'},'Location','north')
    title('Am activation')
subplot(3,1,3)
    plot(alinfo)
    hold on
    plot(100.*Al_musc_all(muscnum,:)')


numinlinks = length(nsys.neuron_objects(muscnum).inlinks);
neuron = nsys.neuron_objects(muscnum);
totdata = zeros(500,1);

% figure('Position',[-(scW-10)/2 130 scW/2 985],'name','inwaves')
% maxSyn = 1000*max([max(.001*(bigH')-.06,[],'all') max(sdata(1).data(dstart:end-10,:),[],'all')]);
% minSyn = 1000*min([min(.001*(bigH')-.06,[],'all') min(sdata(1).data(dstart:end-10,:),[],'all')]);
% for ii = 1:numinlinks
%     neurStr = neuron.inlinks(ii).origin_ID;
%     synNum = str2double(neurStr(regexp(neurStr,'\d')))-38;
%     ndataind = find(contains(sdata(1).data_cols,neurStr(1:end-3)));
%     simdata = sdata(1).data(dstart:end-10,ndataind);
%     calcdata = (bigH(synNum,:)')-60;
%         %nsimdata = 1000.*(interp1(1:length(nsimdata),nsimdata,linspace(1,length(nsimdata),length(mndata)))'+.06);
%         simdata = interp1(1:length(simdata),simdata,linspace(1,length(simdata),length(calcdata)))';
%     %totdata = totdata+(.001*(relW(muscnum,synNum).*bigH(synNum,:)')-.06);
%     subplot(numinlinks,1,ii)
%     plot(calcdata)
%     ylim([minSyn maxSyn])
%     hold on
%     plot(1000*simdata)
%     legend({'Calculated Wave';'Simulation Results'},'Location','north')
%     title([neurStr,'; Synergy ',num2str(synNum)])
% end
% subplot(numinlinks+1,1,ii+1)
% totdata = .001*(relW(muscnum,:)*bigH)-.06;
% plot(totdata)
% hold on
% plot(interp1(1:length(mnsim_results),mnsim_results,linspace(1,length(mnsim_results),length(totdata))))
% %ylim([yMin yMax])
% legend({'Calculated Wave';'Simulation Results'})
% title('Total Calculated Input Voltage')