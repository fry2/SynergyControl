function filetext = inject_joint_waveforms(filetext,equations,end_time)
    % Injects joint waveform equations into joint motor stimuli and overwrites the file. 
    % Input: filepath: simulation .asim file path
    % Input: equations: cell array of sum of sines equations for hip, knee, and ankle
    % Input: end_time: end time of simulation
    
    hip_ind = find(contains(filetext,'<Name>Hipwalking'))+10;
    knee_ind = find(contains(filetext,'<Name>Kneewalking'))+10;
    ankle_ind = find(contains(filetext,'<Name>xAnklewalking'))+10;
    equation_inds = [hip_ind,knee_ind,ankle_ind];
    simtime_ind = find(contains(filetext,'<APIFile/>'))+1;
    
    % Optional: disable joint limits
    limit_enable_ind = contains(filetext,'<EnableLimits>True</EnableLimits>');
    filetext(limit_enable_ind) = deal({'<EnableLimits>False</EnableLimits>'});
    
    filetext{simtime_ind} = number_injector(filetext{simtime_ind},num2str(end_time));
    %project_file{plot_ind} = number_injector(project_file{plot_ind},num2str(end_time));
    
    for i = 1:length(equation_inds)
        old_eq_line = filetext{equation_inds(i)};
        old_endtime_line = filetext{equation_inds(i)-2};
        filetext{equation_inds(i)} = number_injector(old_eq_line,equations{i});
        filetext{equation_inds(i)-2} = number_injector(old_endtime_line,num2str(end_time-.01));
    end
    
    function new_line = number_injector(old_line,new_eq)
        geq = find(old_line=='>');
        leq = find(old_line=='<');
        new_line = [old_line(1:geq(1)),new_eq,old_line(leq(2):end)];
    end
end