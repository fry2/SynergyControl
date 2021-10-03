 function out = wedge(p)
    % wedge operator as defined by Murray Lee and Sastry 1994, p. 41 eq. (2.31), p. 26 (2.4)
    [r,c] = size(p);
%     if ~any(size(p)==3)
%         error('Error wedge.m: Please enter a 3 element vector')
%     end
    if all([r,c]==[3,1]) || all([r,c]==[1,3])
        out = [0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];
    elseif all([r,c]==[6,1]) || all([r,c]==[1,6])
        v = p(1:3);
        w = p(4:6);
        out = [[wedge(w) v];[zeros(1,3) 0]];
    else
        error('Error wedge.m: Size of input must be 3x1 or 4x4')
    end
end