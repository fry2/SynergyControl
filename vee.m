function outB = vee(B)
    % vee operator as defined by MLS 1994 p. 41 (2.30)
    [r,c] = size(B);
    if ~all([r,c]==4)
        error('Error vee: input Matrix not 4x4')
    end
    w1 = B(1:3,1:3);
    v = B(1:3,4);
    outB = [v;w1(3,2);w1(1,3);w1(2,1)];
end