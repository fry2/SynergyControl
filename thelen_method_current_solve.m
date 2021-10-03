function outI = thelen_method_current_solve(A,muscMat)
    solverOptions = optimoptions('fsolve','Display','none');
    outI = zeros(size(muscMat,1),1);
    for ii = 1:size(muscMat,1)
        STmax = muscMat(ii,5); S = muscMat(ii,6); xoff = muscMat(ii,7);
        STcurve = @(I) STmax./(1+exp(S.*(xoff-(I-60)./1000)))-A(ii);
        tempInd = 1;
        for jj = 0:10:20
            if tempInd > 1
               temp = fsolve(STcurve,jj,solverOptions);
               bestTest = STcurve(temp);
               if abs(bestTest) < abs(bestVal)
                   bestI = jj;
               end
            else
               temp = fsolve(STcurve,jj,solverOptions);
               bestVal = STcurve(temp);
               bestI = jj;
            end
            tempInd = tempInd + 1;
        end
        outI(ii) = fsolve(STcurve,bestI,solverOptions);
    end
end