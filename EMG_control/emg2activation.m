function activation = emg2activation(emgInput)

% Coefficients for the recursive filter from Lloyd and Besier
    rowNum = 2; % INPUT
    vals = [-0.033 -0.019 -0.200;...
            -0.091 -0.093 -1.975;...
            +0.265 -0.182 -0.955;...
            -0.097 -0.313 -1.287;...
            +0.006 -0.014 -0.938;...
            +0.015 -0.033 -0.708];

    C1 = vals(rowNum,1);
    C2 = vals(rowNum,2);
    A  = vals(rowNum,3);
    b1 = C1+C2;
    b2 = C1*C2;
    a  = 1+b1+b2;
    u  = zeros(size(emgInput));

    % Apply recursive filter and calculate muscle activation
    for ii = 1:length(emgInput)
        for jj = 1:size(emgInput,2)
            if ii == 1
                u(ii,jj) = a*emgInput(ii,jj);
            elseif ii == 2
                u(ii,jj) = a*emgInput(ii,jj)-b1*u(ii-1,jj);
            else
                u(ii,jj) = a*emgInput(ii,jj)-b1*u(ii-1,jj)-b2*u(ii-2,jj);
            end
        end
    end
    activation = (exp(A.*u)-1)./(exp(A)-1);
end