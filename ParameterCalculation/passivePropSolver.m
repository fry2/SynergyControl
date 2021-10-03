function [ks,kp,stmax] = passivePropSolver(Lr,Fo)
%         beta = 1/(Lmax/Lr-1);
        passPropsolver = @(x) passivePropSolver_func(x,Lr,Fo);
    
    rng default % for reproducibility
    N = 100; % try 100 random start points
    pts = 1000*rand(N,1);
    pts(:,2) = 100*rand(N,1);
    pts(:,3) = 1.05*Fo;
    opts = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');
    %temp_test = @(x,Fo) Fo/x(3)-x(1)/(x(1)+x(2));
    soln = zeros(N,3); % allocate solution
    fval = zeros(N,1); % allocate solution
    parfor k = 1:N
        [soln(k,:),fval(k,:)] = fsolve(passPropsolver,pts(k,:),opts); % find solutions
    end
          
    soln(soln<0) = 0;
    
    sortedsolns = sortrows(soln,3);

    % For the sorted solutions, eliminate those with negative conditional equation outputs
    fval_sort = zeros(length(sortedsolns),1);
    %tempVal = zeros(length(sortedsolns),1);
    for ii = 1:length(sortedsolns)
        fval_sort(ii,1) = passPropsolver(sortedsolns(ii,:));
        %tempVal(ii,1) = temp_test(sortedsolns(ii,:),Fo);
        if fval_sort(ii,1)<0
            sortedsolns(ii,:) = 0;
        end
    end
    
    % Now that only valid solutions exist, pick the one that has the closest STmax to Fmax while still being larger than it

%     for ii = 1:length(sortedsolns)
%         kstemp = sortedsolns(ii,1);
%         kptemp = sortedsolns(ii,2);
%         min_factor(ii,1) = ((kstemp+kptemp)/kstemp)*Fo;
%     end
%     clear kstemp kptemp
%     
%     sortedsolns = sortedsolns(sortedsolns(:,3)>min_factor,:);
    
    if all(sortedsolns==0,'all')
        pts = 1000*rand(N,1);
        pts(:,2) = 100*rand(N,1);
        pts(:,3) = 15*rand(N,1);
        soln = [0 0 0]; % allocate solution
        fTemp = -1;
        %condTemp = -1;
        counter = 1;
        while fTemp < 0 || any(soln<=0) || counter >= 1e4
            kstemp = 1000*rand(1,1);
            kptemp = 100*rand(1,1);
            sttemp = 1.1*Fo;
            soln = fsolve(passPropsolver,[kstemp kptemp sttemp],opts); % find solutions
            fTemp = passPropsolver(soln);
            %condTemp = temp_test(soln,Fo);
            counter = counter + 1;
        end
        ks = kstemp;
        kp = kptemp;
        stmax = sttemp;
    else
        [~,maxKs] = max(sortedsolns(:,1));
        ks = sortedsolns(maxKs,1);
        kp = sortedsolns(maxKs,2);
        stmax = sortedsolns(maxKs,3);
    end   
end