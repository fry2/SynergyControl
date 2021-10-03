function outB = bracket(e1,e2)
    % Bracket operator as defined by MLS 1994 p. 175 (4.26)
    [r1x,c1] = size(e1);
    [r2x,c2] = size(e2);
    if ~all([r1x,r2x]==6) || ~all([c1,c2]==1)
        error('Error bracket: one of the inputs is not 6x1')
    end
    outB = vee(wedge(e1)*wedge(e2)-wedge(e2)*wedge(e1));
end