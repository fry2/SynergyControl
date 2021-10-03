function parameter_injector(fPath,params)
    % For an input project path, update the muscle parameters
    % Input: fPath: Project path to update
    % Input: params: optional cell array with pre-assembled values
    
    %fPath = [pwd,'\Animatlab\SynergyWalking\SynergyWalking20200109_ptest.aproj'];
    
    if ~contains(fPath,'.aproj')
        error('fPath must point to a project file. The input path does not have an .aproj extension.\n')
    end
    
    if ~ischar(fPath)
        fPath = char(fPath);
    end
    
    project_file = importdata(fPath);
    muscle_addresses = contains(project_file,'<Type>LinearHillMuscle</Type>');
    muscle_indices = find(muscle_addresses)-2;
    nummusc = length(muscle_indices);
    
    % If not provided with a param's cell spreadsheet, build one.
    rng(50)
    if nargin == 1
        load([pwd,'\Data\neutral_lengths.mat'],'neutral_lengths');
        lr = neutral_lengths.data(:,1:2);
        load([pwd,'\Data\johnsonMaxForces.mat'],'johnsonMaxForces');
        params = cell(nummusc,4);
        for ii = 1:nummusc
            mInd = muscle_indices(ii);
            lrInd = find(contains(lr(:,1),lower(project_file{mInd}(7:end-7))),1,'first');
            %lrInd = find(contains(strrep(muscle_ranges(:,1),' ',''),lower(project_file{mInd}(10:end-7))),1,'first');
            foInd = find(contains(strrep(johnsonMaxForces(:,1),' ',''),lower(project_file{mInd}(10:end-7))),1,'first');
            Lr = lr{lrInd,2};
            Fo = johnsonMaxForces{foInd,2};
             %rInd = find(contains(project_file(mInd:end),'<RestingLength'),1,'first')-1;
            %fInd = find(contains(project_file(mInd:end),'<MaximumTension'),1,'first')-1;
             %Lr = number_injector(project_file{mInd+rInd},[]);
            %Fo = number_injector(project_file{mInd+fInd},[]);
            params{ii,1} = cell2mat(extractBetween(lower(project_file{mInd}),'<name>','</name>'));
            [ks,kp,stmax,steep,yoff] = parameterSolver(Lr,Fo);
            params{ii,2} = ks;
            params{ii,3} = kp;
            params{ii,4} = stmax; %ST_max values, should be slightly larger than Fmax to ensure Am values are possible.
            params{ii,5} = steep;
            params{ii,6} = yoff;
            [params{:,7}] = deal(-60);
            [params{:,8}] = deal(-40);
            [params{:,9}] = deal(0);
            params{ii,10} = stmax;
            params{ii,11} = 100*(.5*Lr); % 'LT<LowerLimit'
            params{ii,12} = 100*Lr; % '<RestingLength'
            params{ii,13} = 100*(1.5*Lr); % 'LT<UpperLimit'
            params{ii,14} = 100*(.5*Lr); % '<Lwidth'
            params{ii,15} = Fo; % '<MaximumTension'
            [params{:,16}] = deal(0); % 'LT<LowerOutput'
            [params{:,17}] = deal(0); % 'Lt<UpperOutput'
        end
    elseif nargin == 0
        error('Please include a path to an Animatlab project.\n')
    end
    
%     lr = load([pwd,'\Data\neutral_lengths.mat'],'lr');
%     params = import_johnson_data(params,lr);
        
    %These are the parameters represented in the .aproj file. We will iterate through them to update values. The order of these are important and correspond to
    %the column order in the params variable
    parameter_terms = {'<Kse Value';...
                       '<Kpe Value';...
                       '<B Value';...
                       '<C Value';...
                       '<D Value';...
                       'ST<LowerLimit';...
                       'ST<UpperLimit';...
                       'ST<LowerOutput';...
                       'ST<UpperOutput';...
                       'LT<LowerLimit';...
                       '<RestingLength';...
                       'LT<UpperLimit';...
                       '<Lwidth';...
                       '<MaximumTension';...
                       'LT<LowerOutput';...
                       'LT<UpperOutput'};
                    
    for i=1:nummusc
        muscind = muscle_indices(i);
        for j = 1:length(parameter_terms)
            parterm = parameter_terms{j};
            if contains(parameter_terms{j},{'ST';'LT'})
                switch contains(parterm,'LT')
                    case 1
                        curvestr = 'LengthTension';
                    case 0
                        curvestr = 'StimulusTension';
                end
                curvestart = find(contains(project_file(muscind:end),curvestr),1,'first')+muscind-1;
                par_ind = find(contains(project_file(curvestart:end),parterm(3:end)),1,'first')+curvestart-1;
                uselimind = find(contains(project_file(curvestart:end),'UseLimits'),1,'first')+curvestart-1;
                project_file{uselimind} = '<UseLimits>True</UseLimits>';
            else
                par_ind = find(contains(project_file(muscind:end),parameter_terms{j}),1)+muscind-1;
            end
            project_file{par_ind} = number_injector(project_file{par_ind},round(params{i,j+1},2));
        end
    end
    
    [~,projName] = fileparts(fPath);
    carry_on = input(['You are about to overwrite the Animatlab project file (',projName,') you''re using with new parameters.\n'...
                        'This could permanently ruin your project file if there are errors.\n'...
                        'If this is what you want to do, type YES. Otherwise, the program will not overwrite.\n'],'s');
                    
    if strcmp(carry_on,'YES')
        % Record values into a spreadsheet
        spreadsheet = [pwd,'\ParameterCalculation\MuscleParameters_20190513.xlsx'];
        params2write = sortrows(params);
        Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        temp_spr = importdata(spreadsheet);
        spreadsheet_muscles = strrep(temp_spr.textdata(2:end,1),' ','');
        for ii = 1:size(params2write,1)
            sprrow = find(contains(spreadsheet_muscles,params2write{ii,1}(4:end)),1,'first');
            temp_spr.data(sprrow,6:end) = cell2mat(params2write(ii,2:end));
        end
        xlswrite(spreadsheet,temp_spr.data(:,6:end),[Alphabet(7),'2:',Alphabet(7+13),'39']);
        % Update the provided project file with the new parameters
        fileID = fopen(fPath,'w');
        fprintf(fileID,'%s\n',project_file{:});
        fclose(fileID);
    end
    
    function new_line = number_injector(old_line,number)
        quote_inds = find(old_line=='"');
        value = num2str(number);
        if isempty(number) && nargin == 2
            new_line = str2double(cell2mat(extractBetween(old_line,'Actual="','"/>')));
            return
        end
        if contains(old_line,'None')
            actual = value;
        else
            if contains(old_line,'<RestingLength')
                actual = num2str(number/100);
                new_line = [old_line(1:quote_inds(1)),value,old_line(quote_inds(2):quote_inds(3)),'centi',old_line(quote_inds(4):quote_inds(5)),...
                    actual,old_line(quote_inds(6):end)];
                return
            else
                modifier = old_line(quote_inds(3)+1:quote_inds(4)-1);
            end
            switch modifier
                case 'milli'
                    actual = num2str(number/1000);
                case 'centi'
                    actual = num2str(number/100);
            end
        end
        new_line = [old_line(1:quote_inds(1)),value,old_line(quote_inds(2):quote_inds(5)),actual,old_line(quote_inds(6):end)];
    end

    function revised_params = import_johnson_data(params,neutral_lengths)
        % Expand the input parameter array to include information from Johnson 2011 parameters and calculated LT limits/values
        % Input: params: nx6 parameter cell array containing [muscle name, Kse, Kpe, Am, steepness, y offset(?)]
        johnson_excel = [pwd,'\ParameterCalculation\MuscleParameters_20190513.xlsx'];
        C = readcell(johnson_excel);
        [rows,cols] = size(params);
        revised_params = cell(rows,cols+5);
        revised_params(1:rows,1:cols) = params;
        
        %Generate a list of the Johnson data muscle names. These aren't in the same order as the project file
        johnsonforces = cell(length(C),2);
        for q = 2:length(C)
            temp = lower(C{q,1});
            johnsonforces{q,1} = temp(~isspace(temp));
            johnsonforces{q,2} = C{q,6};
        end
        
        clear C johnson_excel
        
        for p = 1:size(revised_params,1)
            %Determine which Johnson muscle to load
            johnsonind = find(contains(johnsonforces(2:end,1),revised_params{p,1}(4:end)),1)+1;
            % Set l_rest to the neutral length
            l_rest = neutral_lengths.lr{p,2}*100;
            % Lower length limit of the LT curve
            revised_params{p,cols+1} = .5*l_rest;
            % L rest
            revised_params{p,cols+2} = l_rest;
            % Upper length limit of the LT curve
            revised_params{p,cols+3} = 1.5*l_rest;
            % L width
            revised_params{p,cols+4} = .5*l_rest;
            % Pull in the optimal force information from the spreadsheet
            f_o = johnsonforces{johnsonind,2};
            % Maximum tension
            revised_params{p,cols+5} = f_o;
        end
    end
end