function [outVal,A] = thelen_method_obj_func(fs,muscMat)
   %solverOptions = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');
   solverOptions = optimoptions('lsqnonlin','Display','none');
   for ii = 1:length(muscMat)
        Ks = muscMat(ii,1); Kp = muscMat(ii,2); Lr = muscMat(ii,3); STmax = muscMat(ii,5);
        solve_func = @(As) (Ks/(Ks+Kp)).*(1+(Kp.*Lr./(4.*As)).^2).*As-fs(ii);
        if fs(ii) == 0
            A(ii) = 0;
        else
            %A(ii) = fsolve(solve_func,fs(ii)/2,solverOptions);
            A(ii) = lsqnonlin(solve_func,fs(ii)/2,0,STmax,solverOptions);
        end
        Arel(ii) = A(ii)/STmax;
   end
   outVal = sum(Arel.^2);
end