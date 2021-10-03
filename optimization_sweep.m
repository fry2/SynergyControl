% Purpose: Distribute forces among leg muscles for a given torque input by
% Sweeping 5 different cost functions and outputting the force profiles for each

close all

vars = struct2cell(whos);

if ~any(contains(vars(1,:),'obj'))
    [obj,sim_file,joints,bodies,joint_limits,joint_profile,sdata,passive_tension] = design_synergy;
end

p = gcp('nocreate');
if isempty(p)
    parpool('local',4);
end

%cost_function_to_use = 3;
toplot = 0;
results_cell = cell(10,5);
results_cell(:,1) = {'';
                     'Forces w/Passive';
                     'Forces No Passive';
                     'Torques w/ Passive';
                     'Torques No Passive';
                     'Moment Arm output';
                     'Fval w/ Passive';
                     'Fval No Passive';
                     'Fmt Waveform';
                     'Opt Time'};

results_cell(1,2:6) = {'Pedotti 1';
                       'Pedotti 2';
                       'Pedotti 3';
                       'Pedotti 4';
                       'Seireg 2'};
                 
for i=1:38
    muscle = obj.musc_obj{i};
    Fmax(i,1) = muscle.max_force;
    %delL(:,1) = max(muscle.muscle_length_profile-muscle.RestingLength,0);
end
                   
for cost_function_to_use = 1:5

    switch cost_function_to_use
        case 1
            fun = @(x) sum(x);
        case 2
            fun = @(x) sum(x.^2);
        case 3
            fun = @(x) sum(x./Fmax);
        case 4
            fun = @(x) sum((x./Fmax).^2);
        case 5
            fun = 'minwork';
    end

    %%
    [forces,tau2,moment_output,fval,telapsed,indiv_torques] = obj.optimize_forces(fun,toplot);
    results_cell{2,cost_function_to_use+1} = forces;
    results_cell{3,cost_function_to_use+1} = tau2;
    results_cell{4,cost_function_to_use+1} = moment_output;
    results_cell{5,cost_function_to_use+1} = fval;
    results_cell{6,cost_function_to_use+1} = telapsed;
   
    
%     for i = 1:5
%         subplot(5,1,i)
%         plot(xnew,results_cell{2,i+1},'LineWidth',1.5)
%         ylabel('Force (N)','FontSize',14)
%         a = get(gca,'XTickLabel');
%         b = get(gca,'YTickLabel');
%         set(gca,'YTickLabel',b,'fontsize',16)
%         if i==1 || i==5
%             switch i
%                 case 1
%                     title('Force Decomposition for Five Optimization Cost Functions','FontSize',20)
%                     set(gca,'XTickLabel',[])
%                 case 5
%                     xlabel('Stride (%)','FontSize',18)
%                     set(gca,'XTickLabel',a,'fontsize',16)
%             end
%         else
%             set(gca,'XTickLabel',[])
%         end
%     end
end
