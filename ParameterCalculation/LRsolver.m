function outData = LRsolver(maxminsim)
    % For an input simulation file where the joints cover a max-min motion (any kind will work), output LT parameters
    
    % Input: maxminsim (char or string): file path to a simulation file that has joint motors moving the joints through their full range
    
    % Example maxminsim locations that was used in the past
    %maxminsim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim";
    sdata = processSimData(maxminsim);

    %%
    lendata = sdata(2).data(1:end-20,:);
    muscnames = sdata(2).data_cols;
    lmins = min(lendata)';
    lmaxs = max(lendata)';
    lrs = .5.*(lmaxs+lmins);
    Lwidth_func = @(perc) sqrt(.25/(1-perc))*(lmaxs-lmins);
    lwidths = Lwidth_func(.7);
    data = [muscnames',num2cell(lrs),num2cell(lwidths),num2cell(lmins),num2cell(lmaxs),num2cell(repmat(.7,38,1))];
    outData = struct();
    outData.data = data;
    outData.data_cols = [{'Muscle Names'},{'Lr'},{'Lw'},{'Lmin'},{'Lmax'},{'Perc'}];
end
