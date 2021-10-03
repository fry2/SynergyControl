function [ks,kp,stmax,steep,yoff] = parameterSolver(lr,fo)
    [ks,kp,stmax] = passivePropSolver(lr,fo);
    [steep,yoff] = stPropSolver(stmax,fo);
end