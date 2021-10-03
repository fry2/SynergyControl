function isbroken = temp_sslentest(obj)

    st_curve = @(Fmax,steepness,xoff,V,yoff) (Fmax./(1+exp(steepness*(xoff-V))))+yoff;
    isbroken = zeros(38,3);

    for mnum = 1:38
        muscle = obj.musc_obj{mnum};
        steepness = muscle.steepness;
        yoff = muscle.y_off;
        xoff = muscle.x_off;
        stmax = muscle.ST_max;
        lr = muscle.RestingLength;
        kp = muscle.Kpe;
        ml = muscle.muscle_length_profile;
        Am = st_curve(stmax,steepness,xoff,-.06,yoff);
        a = Am/(lr*kp);
        isbroken(mnum,1) = ml(end)/lr;
        isbroken(mnum,2) = 1+(1-sqrt(1+16*a^2))/(8*a);
        isbroken(mnum,3) = mean(ml)/lr < (1+(1-sqrt(1+16*a^2))/(8*a));
        if isbroken(mnum,3) == 1
            keyboard
        end
    end

end