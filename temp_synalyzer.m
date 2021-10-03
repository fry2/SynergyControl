clear objmnames
close all

%synData = processSimData([pwd,'\Animatlab\SynergyWalking\SynergyWalking20200109_Standalone.asim'],0);
synData =sdata;
mn = synData(3).data;
tn = synData(2).data;

for ii=1:38
    objmnames{ii,1} = obj.musc_obj{ii}.muscle_name;
end

for ii = 1:length(synData(2).data_cols)
    temp = find(contains(objmnames,synData(2).data_cols(ii)));
    muscles(ii,1) = obj.musc_obj{temp};
    goalforces(:,ii) = forces(:,temp);
end

mnum = 1;

m = length(tn(:,mnum));
n = length(goalforces(:,mnum));

if m ~= n
    goalforcesbig = interp1(1:n,goalforces(:,mnum),linspace(1,n,m));
end

figure
    subplot(3,1,1)
        plot(tn(:,mnum))
        title('Simulation Results')
        xlim([0 size(tn(:,mnum),1)])
        legend(muscles(mnum).muscle_name,'Interpreter','None')
    subplot(3,1,2)
        plot(goalforces(:,mnum))
        title('Desired Forces')
        legend(muscles(mnum).muscle_name,'Interpreter','None')
    subplot(3,1,3)
        plot(tn(:,mnum)-goalforcesbig')
        title('Sim-Desired')
        legend(muscles(mnum).muscle_name,'Interpreter','None')