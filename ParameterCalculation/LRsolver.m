function outData = LRsolver(inSim,rSpread)
    % For an input simulation file where the joints cover a max-min motion (any kind will work), output LT parameters
    
    % Input: maxminsim (char or string): file path to a simulation file that has joint motors moving the joints through their full range
    
    % Example maxminsim locations that was used in the past
    %maxminsim = "G:\My Drive\Rat\SynergyControl\Animatlab\SynergyWalking\SynergyControl_Standalone.asim";
    sdata = processSimData(inSim);
    lenInd = find(contains({sdata.name},'KeyMuscleLen'));
    muscnames = sdata(lenInd).data_cols;
    perc = .7; % Percent value at LT bounds
    %rng(1)
    rng(rSpread)
    % NW neutral is [103.58, 116.73, 138.11] in NW coords
    % In [5.1427, 14.504, 21.863] in Animatlab coords [0.089757, 0.25314, 0.38158]
    
    is_maxmin = 0; % Do you want to set Lr to be half of Lmax+Lmin?
    if is_maxmin
        lendata = sdata(lenInd).data(1:end-20,:);
        lmins = min(lendata)';
        lmaxs = max(lendata)';
        lrs = .5.*(lmaxs+lmins);
        Lwidth_func = @(perc) sqrt(.25/(1-perc))*(lmaxs-lmins);
        lwidths = Lwidth_func(perc);
        data = [muscnames',num2cell(lrs),num2cell(lwidths),num2cell(lmins),num2cell(lmaxs),num2cell(repmat(perc,38,1))];
    else
        lendata = sdata(lenInd).data(1000:end-20,:);
        %rSpread = .2;
        %r = (1-rSpread) + ((1+rSpread)-(1-rSpread)).*rand(1,38);
        r = .7 + (1-.7).*rand(1,38);
        %r = 1 + rSpread.*rand(1,38);
        lrs = (r.*mean(lendata))';
        lmins = .5.*lrs;
        lmaxs = 1.5.*lrs;
        Lwidth_func = @(perc1,Lr) Lr.*sqrt(.25./(1-perc1)); % Based on percentage value when L = 1.5Lr (ex. perc = .7)
        lwidths = Lwidth_func(perc,lrs);
        data = [muscnames',num2cell(lrs),num2cell(lwidths),num2cell(lmins),num2cell(lmaxs),num2cell(repmat(perc,38,1))];
    end
        
    %%
    outData = struct();
    outData.data = data;
    outData.data_cols = [{'Muscle Names'},{'Lr'},{'Lw'},{'Lmin'},{'Lmax'},{'Perc'}];
end
