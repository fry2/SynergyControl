function [steep,yoff] = stPropSolver(stmax,Fo)  
    % For an input Stmax and Fo, output Steepness and Yoffset
    
    % Input: stmax (double): maximum value of the activation function (some mutliple of Fo)
    % Input: Fo (double): maximum tension force in Newtons
    % Output: steep (double): steepness of the ST curve
    % Output: yoff (double): yoffset of the ST curve
    
    % Optimization cost function. Defines two equations for an ST curve relationship
    stsolver = @(x) stPropSolver_func(x);

    pattOpts = optimoptions('patternsearch','MaxTime',1*60,'SearchFcn','MADSPositiveBasis2N','UseCompleteSearch',true,'Display','none');
    steep = patternsearch(stsolver,500,[],[],[],[],0,10e4,[],pattOpts);

    yoff = -stmax/(1+exp(.01*steep));
    
%     V = -.06:.1e-6:-.04;
%     Am = @(V,steep,yoff) stmax./(1+exp(steep.*(-.05-V))) + yoff;
%     bb = Am(V,steep,yoff);figure;plot(V,bb)
end

function outF = stPropSolver_func(x)
    % x(1) = steepness
    % x(2) = y_offset
    % At -40mV, Am should equal STmaxfactor*Fmax
    %F(1) = x(2) + (r*Fo)/(1+exp(-.01*x(1))) - r*Fo;
    %F = r*Fo-x(2)/(1/(1+exp(-.01*x(1)))-.98);
    outF = abs(1/(1+exp(-.01*x))-1/(1+exp(.01*x))-.98);
end