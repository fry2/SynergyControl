function sim_eqn = joint_angle_shifter(shift)
    % For an input shift, modify pre-existing joint angle waveform equations by shifting them in time
    close all
    
    % Equations
    syms t
    joint(1) = 0.2466*sin(0.5468*t-.2808) + 0.046*sin(1.003*t+1.326) + 0.07998*sin(3.77*t+7.8833) + 0.5721*sin(1.879*t+1.5994) + 0.6951*sin(0.2447*t+3.8833);

    joint(2) = 0.1995*sin(0.5822*t-.2136) + 0.08642*sin(0.7588*t +2.4032) + 0.4931*sin(0.209*t -2.2166) + 0.1947*sin(3.768*t +7.6587) + 0.01519*sin(5.651*t+7.9274) + 0.1605*sin(1.875*t + 3.0364) + 0.0206*sin(7.535*t+ 13.454) + 0.003347*sin(2.283*t +2.393);

    joint(3) = 0.0587*sin(0.7246*t + 3.6013) + 0.06009*sin(5.647*t +7.1907) + 0.004938*sin(2.263*t +2.3159) + 0.4367*sin(0.3509*t + 1.3349) + 1.141*sin(0.1553*t -1.5717) + 0.2713*sin(3.768*t + 7.9197) + 0.08073*sin(1.84*t + 3.4415) + 0.02449*sin(7.555*t + 14.2198);
    
    % Anon eq for shifting an array of coefficients
    quick_eqn = @(a,qshift) a(:,1).*qshift+a(:,2);
    sim_eqn = sym(zeros(1,3));

    for j = 1:3
        % Generate a cell array, Num, that contains all coefficients in the equations
        B = regexp(char(vpa(joint(j),5)),'([+|-])|(\d+(\.\d+)?)|(\.\d+)','Match');
        placer = 1;
        numplace = 1;
        aa = [];
        while placer < length(B)
            if ~isempty(B{placer})
                if isnan(str2double(B{placer}))
                    switch B{placer}
                        case '-'
                            Num{numplace,j} = -1*str2double(B{placer+1});
                        case '+'
                            Num{numplace,j} = 1*str2double(B{placer+1});
                    end
                    placer = placer+2;
                else
                    Num{numplace,j}=str2double(B{placer});
                    placer = placer + 1;
                end
            end
            numplace = numplace + 1;
        end
        % With Num in hand, shift the relevant coefficients in time
        aa(:,1) = cell2mat(Num(2:3:length(Num),j));
        aa(:,2) = cell2mat(Num(3:3:length(Num),j));
        bb = quick_eqn(aa,shift);
        bcount = 1;
        for k = 3:3:length(Num)
            Num{k,j} = bb(bcount);
            bcount = bcount + 1;
        end
        % Create animatlab sum of sines equations based on the new coefficients
        simcell = cell(length(Num),2);
        simcell(:,2) = Num(:,j);
        sim_eqn(j) = sum_of_sines_maker(simcell,1);
    end
    % Transcribe equations into text file
    fid = fopen( 'Output_Text/shiftedJointAngles.txt', 'wt' );
    for eqn = 1:3
      fprintf( fid, [char(vpa(sim_eqn(eqn),5)),'\n']);
    end
    fclose(fid);
    
    % Plot original equations and then their shifted forms
    for i = 1:2
        figure
        for j = 1:3
            switch i
                case 1
                    eqn_holder = joint;
                case 2
                    eqn_holder = sim_eqn;
            end
            fplot(eqn_holder(j),[0 10])
            hold on
        end
        ylim([-1.2 .3])
    end
end