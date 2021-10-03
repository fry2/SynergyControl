function [fitresult, gof, sinnum, equations, timeBnds] = sumsinesFit(time, inWave,sinnum)
% Use the curve fitting tool to generate a sum-of-sines representation of the input waveform (r^2>=.97) 
[xData, yData] = prepareCurveData( time, inWave );
opts = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off','TolFun',1e-3);
%opts.Display = 'Off';
temp = [0 0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf 0 0 -Inf];

r2 = 0;
if nargin < 3
    sinnum = 4;
end

    while r2 < .97 && sinnum < 9
        lower = temp(1:sinnum*3);
        typeStr = ['sin',num2str(sinnum)];
        ft = fittype(typeStr);
        opts.Lower = lower;
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
        r2 = gof.rsquare;
        sinnum = sinnum+1;
    end
    
    coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
    equations{1,1} = sum_of_sines_maker(coeffs,1);
    equations{1,2} = sum_of_sines_maker(coeffs,0);
    timeBnds = [time(1), time(end)];
end


