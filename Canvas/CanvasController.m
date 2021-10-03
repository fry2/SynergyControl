classdef CanvasController < handle
    
    properties
        model
        view
        fig
        
        graphic_objects = struct()
        
        listeners = struct()
        
        selection = gobjects(0)
        
        has_current_file = false
        current_file = ''
        saved_to_file = false
        
        saved_to_workspace = false
        
        link_drawing = 'off';
        selection_box = 'off';
        selpos
    end
    
    methods % constructor and destructor
        
        function obj  = CanvasController(model, view)
            
            obj.model = model;
            obj.view = view;
            obj.view.Units = 'pixels';
            if strcmp(view.Parent.Type,'uipanel')
                obj.fig = view.Parent.Parent;
            else
                obj.fig = view.Parent;    
            end
            
            
            obj.listeners.viewDeleted = event.listener(obj.view, ...
                'viewDeleted', @(~,~) delete(obj));
            obj.listeners.modelChanged = event.listener(obj.model, ...
                'modelChanged', @(~,~) obj.modelChanged());
            obj.listeners.graphicItemCreated = event.listener(obj.view, ...
                'graphicItemCreated', ...
                @(~,data) obj.addItemCallback(data.Type, data.Index));
            obj.listeners.canvasSizeChanged = event.listener(obj.view, ...
                'canvasSizeChanged', ...
                @(~,~) obj.setupContextMenu('normal'));
            
            obj.setupMenus();
            obj.setupToolbar();
            obj.setupContextMenu('normal');
            obj.refreshContextMenus();
            
            obj.fig.WindowButtonMotionFcn = @(~,~) obj.defaultMotion();
            obj.fig.WindowButtonDownFcn = @(~,~) obj.defaultClick();
            obj.fig.WindowButtonUpFcn = @(~,~) obj.defaultUnClick();
            obj.fig.WindowKeyPressFcn = @(~,data) obj.defaultKeyPress(data);
        end
        
        function delete(obj)
            
            structfun(@(x) delete(x), obj.listeners);
            structfun(@(x) delete(x), obj.graphic_objects );
            if isvalid(obj.view) && strcmpi(obj.view.BeingDeleted, 'off')
                delete( [obj.view.graphic_objects.Neurons.UIContextMenu] );
                delete( [obj.view.graphic_objects.Obstacles.UIContextMenu] );
            end
            obj.fig.WindowButtonMotionFcn = '';
            obj.fig.WindowButtonDownFcn = '';
            obj.fig.WindowButtonUpFcn = '';
            obj.fig.WindowKeyPressFcn = '';
            
        end
        
    end
    
    methods (Access = private)
        %% setupMenus
        function setupMenus(obj)
            obj.setupFileMenu;
            obj.setupEditMenu;
            obj.setupImportMenu;
            obj.setupExportMenu;
        end
        %% setupFileMenu
        function setupFileMenu(obj)
            obj.graphic_objects.fileMenu = uimenu('Parent', obj.fig, ...
                'Label', 'File');
            uimenu('Parent', obj.graphic_objects.fileMenu,...
                'Label', 'New', ...
                'Accelerator', 'n',...
                'Callback', @(~,~) obj.mNewCanvas());
            uimenu('Parent', obj.graphic_objects.fileMenu,...
                'Label', 'Open',...
                'Accelerator', 'o',...
                'Separator', 'on', ...
                'Callback', @(~,~) obj.mOpen());
            uimenu('Parent', obj.graphic_objects.fileMenu,...
                'Label', 'Save', ...
                'Accelerator', 's',...
                'Callback', @(~,~) obj.mSave);
            uimenu('Parent', obj.graphic_objects.fileMenu,...
                'Label', 'Save As...', ...
                'Callback', @(~,~) obj.mSaveAs);
            uimenu('Parent', obj.graphic_objects.fileMenu,...
                'Label', 'Quit',...
                'Accelerator', 'q', ...
                'Separator', 'on',...
                'Callback', @(~, ~) close(m.fig));
        end
        %% setupEditMenu
        function setupEditMenu(obj)
            
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            obj.graphic_objects.editMenu = uimenu('Parent', obj.fig, ...
                'Label', 'Edit');
            uimenu('Parent', obj.graphic_objects.editMenu,...
                'Label', 'Delete Neuron',...
                'Enable', 'off', ...
                'Callback', @(~,~) obj.deleteSelection());
            uimenu('Parent', obj.graphic_objects.editMenu,...
                'Label', 'Add Neuron',...
                'Separator', 'on',...
                'Callback', @(~,~) obj.model.addItem('n',...
                [randi(bounds(1)-1) ...
                 randi(bounds(2)-1)]));
            uimenu('Parent', obj.graphic_objects.editMenu,...
                'Label', 'Add Adapter',...
                'Callback', @(~,~) obj.model.addItem('a',...
                [randi(bounds(1)-1) ...
                 randi(bounds(2)-1)]));
        end
        %% setupImportMenu
        function setupImportMenu(obj)
            obj.graphic_objects.importMenu = uimenu('Parent', obj.fig, ...
                'Label', 'Import');
            uimenu('Parent', obj.graphic_objects.importMenu,...
                'Label', 'From Workspace', ...
                'Callback', @(~,~) obj.mImportFromWorkspace());
            uimenu('Parent', obj.graphic_objects.importMenu,...
                'Label', 'From File', ...
                'Callback', @(~,~) obj.mImportSubunit());
%             for k=1:2
%                 if k==1
%                     sep = 'on';
%                 else
%                     sep = 'off';
%                 end
%                 uimenu('Parent', obj.graphic_objects.importMenu,...
%                     'Label', ['Scenario ' num2str(k)], ...
%                     'Callback', @(~,~) obj.importScenario(k), ...
%                     'Separator', sep);
%             end
        end
        %% setupExportMenu
        function setupExportMenu(obj)
            obj.graphic_objects.exportMenu = uimenu('Parent', obj.fig, ...
                'Label', 'Export');
            uimenu('Parent', obj.graphic_objects.exportMenu,...
                'Label', 'To Workspace',...
                'Callback', @(~,~) obj.mSaveToWorkspace());
            uimenu('Parent', obj.graphic_objects.exportMenu,...
                'Label', 'To File',...
                'Callback', @(~,~) obj.mSaveSubunit());
        end
        %% setupToolbar
        function setupToolbar(obj)
            obj.graphic_objects.tbar = uitoolbar(obj.fig);
            im = im2double(imread('file_new.png',...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'New Scenario',...
                'ClickedCallback', @(~, ~) obj.mNewCanvas());
            im = im2double(imread('file_open.png', ...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'Open Scenario',...
                'Separator', 'on',...
                'ClickedCallback', @(~, ~) obj.mOpen);
            im = im2double(imread('file_save.png', ...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'Save Scenario',...
                'ClickedCallback', @(~, ~) obj.mSave);
%             im = im2double(imread('Import_16.png',...
%                 'BackGroundColor', [0.94 0.94 0.94]));
%             uipushtool('Parent', obj.graphic_objects.tbar,...
%                 'CData', im, ...
%                 'ToolTipString', 'Import Current Scenario',...
%                 'Separator', 'on',...
%                 'ClickedCallback', @(~, ~) obj.mImportFromWorkspace());
%             im = im2double(imread('Export_16.png',...
%                 'BackGroundColor', [0.94 0.94 0.94]));
%             uipushtool('Parent', obj.graphic_objects.tbar,...
%                 'CData', im, ...
%                 'ToolTipString', 'Save to Workspace and Export as Current Scenario',...
%                 'Separator', 'off',...
%                 'ClickedCallback', @(~, ~) obj.mSaveToWorkspace());
            im = im2double(imread('Import_16.png',...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'Import Subunit',...
                'Separator', 'on',...
                'ClickedCallback', @(~, ~) obj.mImportSubunit());
            im = im2double(imread('Export_16.png',...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'Export Canvas as Subunit',...
                'Separator', 'off',...
                'ClickedCallback', @(~, ~) obj.mExportSubunit());
            im = im2double(imread('animatlab_icon_small.png',...
                'BackGroundColor', [0.94 0.94 0.94]));
            uipushtool('Parent', obj.graphic_objects.tbar,...
                'CData', im, ...
                'ToolTipString', 'Build Animatlab Project',...
                'Separator', 'on',...
                'ClickedCallback', @(~, ~) obj.mCreateAnimatlabProject());
%             im = im2double(imread('SimulinkModel_16.png',...
%                 'BackGroundColor', [0.94 0.94 0.94]));
%             uipushtool('Parent', obj.graphic_objects.tbar,...
%                 'CData', im, ...
%                 'ToolTipString', 'Open Model',...
%                 'Separator', 'on',...
%                 'ClickedCallback', @(~, ~) obj.mOpenModel());
        end
        %% setupContextMenu: Right Clicking on Figure
        %What happens when you right click somewhere on the model. Create a specified object at this location
        function setupContextMenu(obj,contexttype)
                if strcmp(obj.fig.Children(1).Type,'uicontextmenu')
                    delete(obj.graphic_objects.cm);
                    delete(obj.view.UIContextMenu);
    %                 obj.graphic_objects.cm = obj.fig.Children(1);
                end
                obj.graphic_objects.cm = uicontextmenu(obj.fig);
            if strcmp(contexttype,'normal')
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Label', 'Add Neuron',...
                        'Callback', @(~,~) obj.addItemHere('n'));
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Separator','on',...
                        'Label', 'Clear Selections',...
                        'Callback', @(~,~) obj.clearSelections());
                obj.view.UIContextMenu = obj.graphic_objects.cm;
            elseif strcmp(contexttype,'groupselected')
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Label', 'Add Together',...
                        'Callback', @(~,~) obj.addNetworkAddition());
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Label', 'Subtraction Network',...
                        'Callback', @(~,~) obj.addNetworkSubtraction());
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Label', 'Division Network',...
                        'Callback', @(~,~) obj.addNetworkDivision());
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Label', 'Multiplication Network',...
                        'Callback', @(~,~) obj.addNetworkMultiplication());
                uimenu('Parent', obj.graphic_objects.cm,...
                        'Separator','on',...
                        'Label', 'Delete Selected',...
                        'Callback', @(~,~) obj.deleteSelection());
                obj.view.UIContextMenu = obj.graphic_objects.cm;
            end
        end
        %% addItemCallback: Right Clicking a Valid Item
        function addItemCallback(obj, type, model_index)
            
            cmenu = uicontextmenu(obj.fig);
            if type == 'n'
                label1 = 'Add Link';
                label2 = 'Add Stimulus';
                label3 = 'Delete Neuron';
                vect = obj.view.graphic_objects.Neurons;
            end            
            uimenu(cmenu, ...
                'Label', label1, ...
                'Callback', @(~,~) obj.addLink());
            uimenu(cmenu, ...
                'Label', label2, ...
                'Callback', @(~,~) obj.addStimulus());
            uimenu(cmenu, ...
                'Label', label3, ...
                'Callback', @(~,~) obj.deleteSelection(),...
                'Separator','on');
            vect(model_index).UIContextMenu = cmenu;
            
        end
        %% refreshContextMenus
        function refreshContextMenus(obj)
            
            for k=1:obj.model.num_neurons
                obj.addItemCallback('n', k);
            end
            
        end
    end
    
    methods (Access = private)
        %% mNewCanvas
        function mNewCanvas(obj)
            
            answer = inputdlg({'Number of Neurons:',...
                'Number of Obstacles:'}, ...
                'New', 1, {'1', '0'});
            
            if ~isempty(answer)
                bounds(1) = obj.view.graphic_objects.axes.XLim(2);
                bounds(2) = obj.view.graphic_objects.axes.YLim(2);
                numNeurons = str2double(answer{1});
                numObstacles = str2double(answer{2});
                
                obj.model.newCanvas(numNeurons,bounds);
            end
            
        end
        %% mSave: Save the entire Canvas
        function mSave(obj)
            
            if ~obj.has_current_file
                obj.mSaveAs();
            else
                obj.saveFile(obj.current_file);
                obj.saved_to_file  = true;
                obj.updateFigureName();
            end
            
        end
        %% mSaveAs: SaveAs the entire Canvas
        function mSaveAs(obj)
            
            file_dir = fileparts(mfilename('fullpath'));
            file_dir = fullfile(file_dir,...
                CanvasConstants.DATA_RELATIVE_PATH, 'Custom');
            [filename, pathname] = uiputfile(file_dir, ...
                'Save Scenario Configuration');
            if ~isnumeric(filename) && ~isnumeric(pathname)
                obj.saveFile(fullfile(pathname, filename));
                obj.has_current_file = true;
%                 obj.current_file = filename;
            if ~strcmp(filename(end-1:end),'.m')
                filename = [filename(1:end-7), '.m'];
            end
                obj.current_file = fullfile(pathname, filename);
                obj.saved_to_file = true;
                obj.updateFigureName();
            end
            
        end
        %% mOpen: Open an existing Canvas
        function mOpen(obj)
            
            file_dir = fileparts(mfilename('fullpath'));
            file_dir = fullfile(file_dir,...
                CanvasConstants.DATA_RELATIVE_PATH, 'Custom');
            [filename, pathname] = uigetfile('*.m', ...
                'Open Scenario Configuration', file_dir);
            if ~isnumeric(filename) && ~isnumeric(pathname)
                obj.openFile(fullfile(pathname,filename));
                obj.has_current_file = true;
                obj.current_file = fullfile(pathname,filename);
                obj.saved_to_file = true;
                obj.updateFigureName();
            end
            
        end
        %% mImportSubunit: Import a Subunit from a file
        function mImportSubunit(obj)
            file_dir = fileparts(mfilename('fullpath'));
            file_dir = fullfile(file_dir,...
                CanvasConstants.DATA_RELATIVE_PATH,'\user_networks');
            [filename, pathname] = uigetfile('*.m','Open Scenario Configuration', file_dir);
            obj.openFile([pathname,filename]);
            obj.saved_to_workspace = true;
            obj.saved_to_file = false;
            obj.has_current_file = false;
            obj.current_file = '';
            obj.updateFigureName();
        end
        %% mImportFromWorkspace: Import information from the Workspace
        function mImportFromWorkspace(obj)
            
            file_dir = fileparts(mfilename('fullpath'));
            filename = fullfile(file_dir,...
                CanvasConstants.DATA_RELATIVE_PATH, ...
                CanvasConstants.CURRENT_CANVAS_SETUP);
            obj.openFile(filename);
            obj.saved_to_workspace = true;
            obj.saved_to_file = false;
            obj.has_current_file = false;
            obj.current_file = '';
            obj.updateFigureName();

            
        end
        %% mExportSubunit: Save the current canvas as a subunit to a file
        function mExportSubunit(obj)
            
            file_dir = fileparts(mfilename('fullpath'));
            file_dir = fullfile(file_dir,CanvasConstants.DATA_RELATIVE_PATH,'\user_networks');
            [filename, pathname] = uiputfile(file_dir,'Save Subunit');
            obj.saveSubunit([pathname,filename]);
            
            S = obj.model.getData();
            
            assignin('base','Neuron_Objects', S.Neuron_Objects);
            assignin('base','Neuron_Objects', S.Link_Objects);
            
            obj.saved_to_workspace = true;
            obj.updateFigureName();
            
        end
        %% mSaveToWorkspace: save the existing canvas to the workspace
        function mSaveToWorkspace(obj)
            
            filename = fileparts(mfilename('fullpath'));
            filename = fullfile(filename, ...
                CanvasConstants.DATA_RELATIVE_PATH, ...
                CanvasConstants.CURRENT_CANVAS_SETUP);
            obj.saveFile(filename);
            
            S = obj.model.getData();
            
            % Export in base workspace
            assignin('base','ObstaclesPositions', S.ObstaclesPositions);
            assignin('base','ObstaclesVertices', S.ObstaclesVertices);
            assignin('base','ObstaclesEdges', S.ObstaclesEdges);
            assignin('base','NeuronsPositions', S.NeuronsPositions);
            assignin('base','Neuron_Objects', S.Neuron_Objects);
            
            obj.saved_to_workspace = true;
            obj.updateFigureName();
            
        end
        %% mCreateAnimatlabProject: generate an animatlab project file from canvas
        function mCreateAnimatlabProject(obj)
            obj.model.create_animatlab_project;
        end
        %% saveFile: save SUBUNIT file in \user_networks
        function saveSubunit(obj, filename)
            
            S = obj.model.getData();
            
            % Open new file for writing
            [pathstr, name, ~] = fileparts(filename);
            %Rename the file as a .m since the default save type is .mldatx
            fid = fopen([fullfile(pathstr, name) '.m'],'wt');
            
            % Set file header (to write that file is autogenerated and print the
            % date)
            fprintf(fid,'%% File Autogenerated by function %s\n',mfilename);
            fprintf(fid,'%% Date %s\n',datestr(now,'dddd dd mmmm yyyy HH:MM:SS'));
            
            positions = obj.model.neurons_positions;
            positionsUL = positions-CanvasConstants.NEURON_size/2;
            positionsBR = positions+CanvasConstants.NEURON_size/2;
            a = max(positionsBR(:,1))-min(positionsUL(:,1));
            b = max(positionsBR(:,2))-min(positionsUL(:,2));
            
            fprintf(fid,'\nsize = [%d , %d];\n\n',a,b);
            
            S_obj = fieldnames(S);
            full_text = {};
            for k = 1:numel(S_obj)
                if k>1
                    obj_header = ['\n\n',S_obj{k},' = cell2struct(...\n'];
                else
                    obj_header = [S_obj{k},' = cell2struct(...\n'];
                end
                fieldNames = fieldnames(S.(S_obj{k}));
                field_text = {'{...\n'};
                content = {obj_header;'[...\n';};
                for j = 1:numel(S.(S_obj{k}))
                    content = [content;'{...\n'];
                    for i = 1:numel(fieldNames)
                        fieldcontent = S.(S_obj{k})(j).(fieldNames{i});
                        field_text = [field_text;['\t''',fieldNames{i},''',...\n']];
                        if ischar(fieldcontent)   
                            content = [content;['\t''',S.(S_obj{k})(j).(fieldNames{i}),''';...\n']];
                        elseif isempty(S.(S_obj{k})(j).(fieldNames{i}))
                            content = [content;'\t[];...\n'];
                        elseif iscell(fieldcontent)
                            content = [content;'\t{...\n'];
                            for ii =1:length(fieldcontent)
                                content = [content;'\t''',fieldcontent{ii},''';...\n'];
                            end
                            content = [content;'};...\n'];
                        elseif isstruct(fieldcontent)
                            if strcmp(S_obj{k}(1:6),'Neuron') && strcmp(fieldcontent(1).ID(1:4),'link')
                                content = [content;'[...\n'];
                                for ii = 1:numel(fieldcontent)
                                    linkID_string = fieldcontent(ii).ID;
                                    link_num_ind = regexp(linkID_string,'\d');
                                    link_num = linkID_string(link_num_ind);
%                                     content = [content;'\t',S_obj{1},'(',link_num,');...\n'];
                                    content = [content;'\t',S_obj{1},'(',num2str(ii),');...\n'];
                                end
                                content = [content;'];...\n'];
                            end
                        else
                            content = [content;S.(S_obj{k})(j).(fieldNames{i})];
                        end
                    end
                    content = [content;'},...\n'];
                end
                    content = [content(1:end,:);'],...\n'];
                    field_text = [field_text(1:numel(fieldNames)+1,:);'},...\n';'1);'];
                    full_text = [full_text;content;field_text];
            end
%                     fprintf(fid,obj_header);
                    for k = 1:numel(full_text)
                        print_seg = full_text{k,1};
                        if ischar(print_seg)
                            fprintf(fid,print_seg);
                        else
                            if sum(size(print_seg)) == 2
                                fprintf(fid,'\t%.3f;...\n',print_seg);
                            else
                                fprintf(fid,'\t[%4.0f,%4.0f];...\n',print_seg);
                            end
                        end
                    end
            
            % Close file
            fclose(fid);
            
        end
        %% openFile: open a SUBUNIT file
        function openFile(obj, filename)
            
            clear Neuron_Objects Link_Objects
            
            if sum(filename) == 0
                return
            end
            
            clear(filename);
            run(filename);

            numNeurons = obj.model.num_neurons;
            numAddNeurons = numel(Neuron_Objects);
            
            numLinks = obj.model.num_links;
            numAddLinks = numel(Link_Objects);
            
            neuronfieldnames = fieldnames(Neuron_Objects(1));
            skipfields = {'name';'ID';'outlinks';'inlinks';'outlink_IDs';'inlink_IDs';'location';'nsize';'ca_act_ID';'ca_deact_ID'};
            
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            for i = (1+numNeurons):(numNeurons+numAddNeurons)
                pos = Neuron_Objects(i-numNeurons).location;
                obj.model.addItem('n',pos, bounds);
                %obj.model.neuron_objects(i) = Neuron_Objects(i-numNeurons);
            end
            
            for i=1:numAddNeurons
                index = i+numNeurons;
                for j = 1:length(neuronfieldnames)
                    if sum(contains(skipfields,neuronfieldnames(j))) == 0
                        obj.model.neuron_objects(index).(neuronfieldnames{j}) = Neuron_Objects(i).(neuronfieldnames{j});
                    end
                end
            end
            
            neuronlist = cell(size(Neuron_Objects),1);
            neuronlist = {Neuron_Objects(:).name};
%             neuronlist = {obj.model.neuron_objects(:).name};
            for i = (1+numLinks):(numLinks+numAddLinks)
                link = Link_Objects(i-numLinks);
                startpos = link.start;
                enndpos = link.end;
                start_ind = find(contains(neuronlist,link.origin_ID(1:end-3)),1)+numNeurons;
                end_ind = find(contains(neuronlist,link.destination_ID(1:end-3)),1)+numNeurons;
                obj.model.addLink(start_ind,end_ind,startpos,enndpos,link.synaptictype);
            end
            
            obj.model.update_link_cdata();
        end
        %% addItemHere: add an item at specific point on canvas
        function addItemHere(obj, type)
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            obj.model.addItem(type, obj.view.CurrentPoint(1,1:2),bounds);
            
        end
        %% addStimulus: add a stimulus to a neuron
        function addStimulus(obj)
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            selPos = obj.selection.Position(1:2)+CanvasConstants.NEURON_size/2;
            obj.model.addItem('stimulus',selPos,bounds)
        end
        %% addLink: add a link between two items
        function addLink(obj)
            if length(obj.selection) == 1
                obj.link_drawing = 'on';
            end
        end
        %% addNetworkAddition: create an addition network from selected objects
        function addNetworkAddition(obj)
            sel = obj.selection;
            if length(sel) == 2
                bounds(1) = obj.view.graphic_objects.axes.XLim(2);
                bounds(2) = obj.view.graphic_objects.axes.YLim(2);
                model_n_pos = obj.model.neurons_positions;
                
                numSelected = length(sel);
                posmat = reshape([sel.Position],[4,size(sel,2)])';
                posmat = posmat(:,1:2)+CanvasConstants.NEURON_size/2;
                midpoint = sum(posmat)./2;
                mid_dist = norm(posmat(1,:)-midpoint);
                offset = -atan(diff(posmat(:,2))/diff(posmat(:,1)));
                output_pos = [midpoint(1)-mid_dist*sin(offset) midpoint(2)-mid_dist*cos(offset)];
                if min(sum(abs(model_n_pos-output_pos),2)) < 5
                    random_change = .5*randi(6,1);
                    output_pos = [midpoint(1)-random_change*mid_dist*sin(offset) midpoint(2)-random_change*mid_dist*cos(offset)];
                end
                obj.model.addItem('n',output_pos, bounds);
                
                model_n_pos = obj.model.neurons_positions;
                sel_mindices = zeros(length(posmat),1);
                for jj = 1:length(posmat)
                    [~,sel_mindices(jj)] = min(sum(abs(model_n_pos-posmat(jj,:)),2));
                end

                [~,output_index] = min(sum(abs(model_n_pos-output_pos),2));
                
                for ii = 1:numSelected
                    obj.model.addLink(sel_mindices(ii),output_index,model_n_pos(sel_mindices(ii),:),model_n_pos(output_index,:),'SignalTransmission1');
                end
            end
            obj.clearSelections();
        end
        %% addNetworkSubtraction: create a subtraction network from selected objects
        function addNetworkSubtraction(obj)
            sel = obj.selection;
            numSelected = length(sel);
            bounds = [obj.view.graphic_objects.axes.XLim(2) obj.view.graphic_objects.axes.YLim(2)];
%             bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            if numSelected > 2
                obj.clearSelections();
                return
            else
                sel = obj.selection;
            end
            list = cell(numSelected,1);
            index_list = zeros(numSelected,1);
            for ii = 1:numSelected
                selpos = sel(ii).Position(1:2)+CanvasConstants.NEURON_size/2;
                [~,model_ind] = max(sum(ismember(obj.model.neurons_positions,selpos(1:2)),2));
                list(ii) = {obj.model.neuron_objects(model_ind).name};
                index_list(ii) = model_ind;
            end
            [indx,~] = listdlg('Name','Subtraction Neuron...',...
                                'ListSize',[200 100],...
                                'PromptString','Which neuron is being subtracted?',...
                                'SelectionMode','single',...
                                'ListString',list);
            not_sel = 3-indx;
            
            posmat = reshape([sel.Position],[4,size(sel,2)])';
            posmat = posmat(:,1:2)+CanvasConstants.NEURON_size/2;
            midpoint = sum(posmat)./2;
            mid_dist = norm(posmat(1,:)-midpoint);
            offset = -atan(diff(posmat(:,2))/diff(posmat(:,1)));
            output_pos = [midpoint(1)-mid_dist*sin(offset) midpoint(2)-mid_dist*cos(offset)];
            obj.model.addItem('n',output_pos, bounds);

            model_n_pos = obj.model.neurons_positions;
%                 sel_mindices = zeros(length(posmat),1);
%                 for jj = 1:length(posmat)
%                     [~,sel_mindices(jj)] = min(abs(sum(model_n_pos-posmat(jj,:),2)));
%                 end

            [~,output_index] = min(sum(abs(model_n_pos-output_pos),2));

            obj.model.addLink(index_list(indx),output_index,model_n_pos(index_list(indx),:),model_n_pos(output_index,:),'SignalModulation2');
            obj.model.addLink(index_list(not_sel),output_index,model_n_pos(index_list(not_sel),:),model_n_pos(output_index,:),'SignalTransmission1');
            obj.clearSelections();
        end
        %% addNetworkDivision: create a division network from selected objects
        function addNetworkDivision(obj)
            sel = obj.selection;
            numSelected = length(sel);
            bounds = [obj.view.graphic_objects.axes.XLim(2) obj.view.graphic_objects.axes.YLim(2)];
%             bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            if numSelected > 2
                obj.clearSelections();
                return
            else
                sel = obj.selection;
            end
            list = cell(numSelected,1);
            index_list = zeros(numSelected,1);
            for ii = 1:numSelected
                selpos = sel(ii).Position(1:2)+CanvasConstants.NEURON_size/2;
                [~,model_ind] = max(sum(ismember(obj.model.neurons_positions,selpos(1:2)),2));
                list(ii) = {obj.model.neuron_objects(model_ind).name};
                index_list(ii) = model_ind;
            end
            [indx,~] = listdlg('Name','Division Neuron...',...
                                'ListSize',[200 100],...
                                'PromptString','Which neuron is the divisor?',...
                                'SelectionMode','single',...
                                'ListString',list);
            not_sel = 3-indx;
            
            posmat = reshape([sel.Position],[4,size(sel,2)])';
            posmat = posmat(:,1:2)+CanvasConstants.NEURON_size/2;
            midpoint = sum(posmat)./2;
            mid_dist = norm(posmat(1,:)-midpoint);
            offset = -atan(diff(posmat(:,2))/diff(posmat(:,1)));
            output_pos = [midpoint(1)-mid_dist*sin(offset) midpoint(2)-mid_dist*cos(offset)];
            obj.model.addItem('n',output_pos, bounds);

            model_n_pos = obj.model.neurons_positions;
%                 sel_mindices = zeros(length(posmat),1);
%                 for jj = 1:length(posmat)
%                     [~,sel_mindices(jj)] = min(abs(sum(model_n_pos-posmat(jj,:),2)));
%                 end

            [~,output_index] = min(sum(abs(model_n_pos-output_pos),2));

            obj.model.addLink(index_list(indx),output_index,model_n_pos(index_list(indx),:),model_n_pos(output_index,:),'SignalModulation3');
            obj.model.addLink(index_list(not_sel),output_index,model_n_pos(index_list(not_sel),:),model_n_pos(output_index,:),'SignalTransmission1');
            obj.clearSelections();
        end
        %% addNetworkMultiplication: create a multiplication network from selected objects
        function addNetworkMultiplication(obj)
            sel = obj.selection;
            numSelected = length(sel);
            bounds = [obj.view.graphic_objects.axes.XLim(2) obj.view.graphic_objects.axes.YLim(2)];
%             bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            if numSelected > 2
                obj.clearSelections();
                return
            else
                sel = obj.selection;
            end
%             list = cell(numSelected,1);
%             index_list = zeros(numSelected,1);
%             for ii = 1:numSelected
%                 selpos = sel(ii).Position(1:2)+CanvasConstants.NEURON_size/2;
%                 [~,model_ind] = max(sum(ismember(obj.model.neurons_positions,selpos(1:2)),2));
%                 list(ii) = {obj.model.neuron_objects(model_ind).name};
%                 index_list(ii) = model_ind;
%             end
%             [indx,~] = listdlg('Name','Division Neuron...',...
%                                 'ListSize',[200 100],...
%                                 'PromptString','Which neuron is the divisor?',...
%                                 'SelectionMode','single',...
%                                 'ListString',list);
%             not_sel = 3-indx;
            
            posmat = reshape([sel.Position],[4,size(sel,2)])';
            posmat = posmat(:,1:2)+CanvasConstants.NEURON_size/2;
            midpoint = sum(posmat)./2;
            mid_dist = norm(posmat(1,:)-midpoint);
            offset = -atan(diff(posmat(:,2))/diff(posmat(:,1)));
            output_pos = [midpoint(1)-mid_dist*sin(offset) midpoint(2)-mid_dist*cos(offset)];
            interneuron_pos = [midpoint(1)-mid_dist*sin(offset+45*(pi/180)) midpoint(2)-mid_dist*cos(offset+45*(pi/180))];
            obj.model.addItem('n',output_pos, bounds);
            obj.model.addItem('n',interneuron_pos, bounds);

            model_n_pos = obj.model.neurons_positions;
            sel_mindices = zeros(length(posmat),1);
            for jj = 1:length(posmat)
                [~,sel_mindices(jj)] = min(sum(abs(model_n_pos-posmat(jj,:)),2));
            end

            [~,output_index] = min(sum(abs(model_n_pos-output_pos),2));
            [~,interneuron_index] = min(sum(abs(abs(model_n_pos-interneuron_pos)),2));

            obj.model.addLink(sel_mindices(1),interneuron_index,model_n_pos(sel_mindices(1),:),model_n_pos(interneuron_index,:),'SignalModulation3');
            obj.model.addLink(interneuron_index,output_index,model_n_pos(interneuron_index,:),model_n_pos(output_index,:),'SignalModulation3');
            obj.model.addLink(sel_mindices(2),output_index,model_n_pos(sel_mindices(2),:),model_n_pos(output_index,:),'SignalTransmission1');
            obj.clearSelections();
        end
        %% importScenario
%         function importScenario(obj, num)
%             
%             file_dir = fileparts(mfilename('fullpath'));
%             filename = fullfile(file_dir,...
%                 CanvasConstants.DATA_RELATIVE_PATH, ...
%                 'Scenarios',...
%                 ['Scenario' num2str(num) '.m']);
%             obj.openFile(filename);
%             obj.saved_to_workspace = true;
%             obj.saved_to_file = false;
%             obj.has_current_file = false;
%             obj.current_file = '';
%             obj.updateFigureName();
%             
%         end
        %% defaultMotion: function called repeatedly as cursor moves across canvas
        function defaultMotion(obj)
            %This is constant, be careful setting a breakpoint.
            %This is really just the hover action
            cp = obj.view.CurrentPoint(1,1:2);
            Npos = obj.model.neurons_positions;
            overN = ~isempty(Npos) && any(all(cp>(Npos-CanvasConstants.NEURON_size/2),2) & ...
                                          all(cp<(Npos+CanvasConstants.NEURON_size/2),2));
              overO = 0;
            if  overO || overN
                obj.fig.Pointer = 'fleur';
            else
                obj.fig.Pointer = 'arrow';
            end
            
            if strcmp(obj.link_drawing,'on')
                sel_pos = obj.selpos(1:2)+CanvasConstants.NEURON_size/2;
                coords = [sel_pos;cp];
                if strcmp(obj.view.graphic_objects.axes.Children(1).Type,'line')
                    delete(obj.view.graphic_objects.axes.Children(1));
                end
                line(coords(:,1),coords(:,2),...
                    'Tag','linkline',...
                    'Parent',obj.view.graphic_objects.axes);                
            end
        end
        %% defaultClick: function called when canvas is clicked
        function defaultClick(obj)
            link_on = 0;
            cp = obj.view.CurrentPoint(1,1:2);
             if isvalid(obj.selection)
                if strcmp(obj.link_drawing,'on')
                    obj.link_drawing = 'off';
                    link_on = 1;
                    Npos = obj.model.neurons_positions;
                    end_neur = [all(cp>(Npos-CanvasConstants.NEURON_size/2),2),all(cp<(Npos+CanvasConstants.NEURON_size/2),2)];
                    end_ind = find(all(end_neur,2));
                    linkpos1 = obj.selection.Position(1:2)+CanvasConstants.NEURON_size/2;
                    start_neur = [all(linkpos1>(Npos-CanvasConstants.NEURON_size/2),2),all(linkpos1<(Npos+CanvasConstants.NEURON_size/2),2)];
                    start_ind = find(all(start_neur,2));
                    if isempty(end_ind) || start_ind==end_ind
                        obj.fig.CurrentObject = [];
                    else
                        obj.fig.CurrentObject = obj.view.graphic_objects.Neurons(1,end_ind);
                    end
                    obj.selection.Selected = 'off';
                end
             end
            
            if sum(contains({obj.view.graphic_objects.Neurons.Selected},'on')) > 1
                %multiple neurons have been selected by a bouding box
                %in order to switch the selection to a neuron outside of that bounded group
                %we need to decide whether the newly selected neuron is in the group or not
                viewNeurPos = cell2mat({obj.view.graphic_objects.Neurons.Position}')+[CanvasConstants.NEURON_size/2 0 0];
                [~,clickedNeur] = min(sum(abs(viewNeurPos(:,1:2)-cp),2));
                selvec = contains({obj.view.graphic_objects.Neurons.Selected},'on');
                if ~selvec(clickedNeur)
                    %the newly clicked neuron is not in the selected group
                    %change the selection to it and deselect the group
                    sel = obj.view.graphic_objects.Neurons(clickedNeur);
                else
                    %the selected neuron is a part of the selected group
                    %keep the group selection the same
                    sel = obj.view.graphic_objects.Neurons(selvec);
                end
            else
                %when not worried about groups, just make the selected object the fig.CurrentObject
                sel = obj.fig.CurrentObject;
                selvec = [];
            end
            
            %include an option that keeps a neuron selected when you select a stimulus tab
            %previously, selecting a tab would deselect the neuron
            if ~isempty(sel) && isempty(selvec) && size(obj.selection,2) < 2
                if ~strcmp(sel.Type,'image') && ~isempty(obj.selection) %"if not selecting the background space and the obj.selection is not empty"
                    if (strcmp(obj.selection.Type,'rectangle') && strcmp(sel.Type,'uitab')) %"if the obj.selection is a neuron and what you just clicked is a tab"
                        sel = obj.selection; %set the selection object back to the neuron, essentially ignoring the tab selection
                    end
                end
            end
            
            if (isa(sel, 'matlab.graphics.primitive.Rectangle') && size(sel,2) == 1)

                obj.setSelection(sel);
             
                sel.Selected = 'on';
                %Find the center of the rectangle
                pos = sel.Position(1:2) + sel.Position(3:4)/2;

                if sel.Tag == 'n'
                    model_n_pos = obj.model.neurons_positions;
                    [~,model_index] = min(sum(abs(model_n_pos-pos),2));
                    obj.selection.UserData{1} = sel.Tag;
                    obj.selection.UserData{2} = model_index;
                    obj.selection.UserData{3} = model_n_pos;
                else
                    model_index = find(all(pos==obj.model.obstacles_positions,2));
                end
                
                obj.view.formSelection(sel.Tag,model_index);
                startCP = obj.view.CurrentPoint(1,1:2);
                if link_on
                    link_on = 0;
                    numAxesChildren = size(obj.view.graphic_objects.axes.Children,1)-1;
                    objtypes = cell(1,numAxesChildren);
                    for i=1:numAxesChildren
                        objtypes{i} = obj.view.graphic_objects.axes.Children(i).Type;
                    end
                    lineind = find(contains(objtypes,'line'));
                    if size(lineind,2) ~= size(obj.model.link_ends,1) + 1
                        link2delete = size(obj.model.link_ends,1);
                        obj.model.deleteItem('l', link2delete);
                    end
                    delete(obj.view.graphic_objects.axes.Children(lineind(1)));
                    if ~all(obj.model.neuron_objects(start_ind).location == linkpos1) || ~all(obj.model.neuron_objects(model_index).location == pos)
                        fprintf('The link end positions don''t match one of the end neuron positions.\n')
                        keyboard
                    end
%                     list = {'Depolarizing IPSP','Hyperpolarizing IPSP'};
                    list = {obj.model.synapse_types.name};
                    [indx,tf] = listdlg('Name','Synpase Types',...
                                        'ListSize',[200 100],...
                                        'PromptString','Select a synapse type',...
                                        'SelectionMode','single',...
                                        'ListString',list);
                    if tf
                        linktype = list{indx};
                        obj.model.addLink(start_ind,model_index,linkpos1,pos,linktype);
                    end
                else
                    obj.fig.WindowButtonMotionFcn = @(~,~) obj.clickedMotion(...
                    sel.Tag, model_index, pos - startCP);
                end

            elseif isa(sel, 'matlab.graphics.primitive.Rectangle') && size(sel,2) > 1
                obj.setSelection(sel);             
                [sel.Selected] = deal('on');
                
                posmat = reshape([sel.Position],[4,size(sel,2)])';
                posmat = posmat(:,1:2)+CanvasConstants.NEURON_size/2;

                [~,selInd] = min(sum(abs(posmat-cp),2));
                pos = posmat(selInd,:);
                
                startCP = obj.view.CurrentPoint(1,1:2);
                
                obj.fig.WindowButtonMotionFcn = @(~,~) obj.clickedMotion(...
                    'group', selInd, pos - startCP);
            elseif isa(sel, 'matlab.graphics.primitive.Line') || isa(sel, 'matlab.graphics.chart.primitive.Quiver')
                model_index = contains({obj.model.link_objects.ID},sel.Tag);
                obj.view.formSelection(sel.Tag,model_index);
            elseif strcmp(obj.fig.SelectionType,'normal')
                obj.emptySelection(link_on);
                obj.fig.WindowButtonMotionFcn = @(~,~) obj.clickedMotionSelBox(cp);
            else
                obj.setupContextMenu('normal');
            end
            
        end
        %% clickedMotionSelBox: creating a group selection
        function clickedMotionSelBox(obj,cp)
            obj.selection_box = 'on';
            bb = 1;
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            
            childData = {};
            for i=1:size(obj.view.graphic_objects.axes.Children,1)
                child = obj.view.graphic_objects.axes.Children(i);
                if isa(child,'matlab.graphics.primitive.Rectangle')
                    if strcmp(child.Tag,'n')
                        dataterm = 'neuron';
                    elseif strcmp(child.UserData,'selectionbox')
                        dataterm = 'selectionbox';
                    end
                else
                    dataterm = 'notarect';
                end
                childData = [childData,dataterm];
            end
%             childData = {obj.view.graphic_objects.axes.Children(1:end-1).UserData};
%             childData = childData(~cellfun(@isempty,childData));
            selInd = find(contains(childData,'selectionbox'));
            
            sp = round(cp);
            hp = obj.view.CurrentPoint(1,1:2);
            hp = obj.model.ConstrainedPosition('selectionbox', hp, bounds);
            
            if all(hp>sp)
                rectpos = [sp(1) sp(2) abs(hp(1)-sp(1)) abs(sp(2)-hp(2))];
            elseif all(hp<sp)
                rectpos = [hp(1) hp(2) abs(hp(1)-sp(1)) abs(sp(2)-hp(2))];
            elseif hp(1) >= sp(1) && hp(2) <= sp(2)
                rectpos = [sp(1) hp(2) abs(hp(1)-sp(1)) abs(sp(2)-hp(2))];
            elseif hp(1) <= sp(1) && hp(2) >= sp(2)
                rectpos = [hp(1) sp(2) abs(hp(1)-sp(1)) abs(sp(2)-hp(2))];
            end
            
            if ~isempty(selInd)
                if ~all(obj.view.graphic_objects.axes.Children(selInd).Position == rectpos)
                     set(obj.view.graphic_objects.axes.Children(selInd),'Position',rectpos);
                end
            else
                rectangle(obj.view.graphic_objects.axes,'Position',rectpos,'UserData','selectionbox','EdgeColor',[0 0 1],'LineWidth',.3);
                drawnow
            end
        end
        %% clickedMotion: action taken when item is clicked and dragged
        function clickedMotion(obj, type, index, delta)
            
            bounds(1) = obj.view.graphic_objects.axes.XLim(2);
            bounds(2) = obj.view.graphic_objects.axes.YLim(2);
            cp = obj.view.CurrentPoint(1,1:2);
            boundpos = bounds-CanvasConstants.NEURON_size/2-[1 1];
            
            if strcmp(type,'group')
                neurpos = obj.model.neurons_positions;
                anchornode = obj.selection(index).Position(1:2)+CanvasConstants.NEURON_size/2;
                anchorrel = neurpos-anchornode;
                nodepos = obj.selection(index).Position(1:2)+CanvasConstants.NEURON_size/2;
                [~,modelInd] = min(sum(abs(neurpos-nodepos),2));
                obj.model.moveItem('n', modelInd, cp+round(delta)+anchorrel(modelInd,:), bounds);
                moved_pos = obj.model.neuron_objects(modelInd).location;
                is_constrained = min(abs(moved_pos-boundpos))<3 || min(abs(moved_pos-[41 41]))<3;
                for i = 1:size(obj.selection,2)
                    if i~=index
                        nodepos = obj.selection(i).Position(1:2)+CanvasConstants.NEURON_size/2;
                        [~,modelInd] = min(sum(abs(neurpos-nodepos),2));
                        mover = cp+round(delta)+anchorrel(modelInd,:);
                        if is_constrained
                            nomove = min(abs([moved_pos-boundpos;moved_pos-[41 41]]))<3;
                            mover(nomove) = obj.model.neurons_positions(modelInd,nomove);
                        end
                        obj.model.moveItem('n', modelInd, mover, bounds);
                    end
                end
            else
                obj.model.moveItem(type, index, cp + delta, bounds);
            end
            
        end
        %% defaultUnClick: callback for when mouse is "unclicked"
        function defaultUnClick(obj)

            if strcmp(obj.selection_box,'on')
                numAxesChildren = size(obj.view.graphic_objects.axes.Children,1);
                childData = cell(1,numAxesChildren);
%                 selCell = {};
                for i=1:numAxesChildren
                    child = obj.view.graphic_objects.axes.Children(i);
                    if isa(child,'matlab.graphics.primitive.Rectangle')
                        if strcmp(child.Tag,'n')
                            dataterm = 'neuron';
                        elseif strcmp(child.UserData,'selectionbox')
                            dataterm = 'selectionbox';
                        end
                    else
                        dataterm = 'notarect';
                    end
                    childData{i} = dataterm;
                end
                selInd = contains(childData,'selectionbox');
                selBoxArea = obj.view.graphic_objects.axes.Children(selInd).Position;
                xbounds = [min(selBoxArea(1),(selBoxArea(1)+selBoxArea(3))) max(selBoxArea(1),(selBoxArea(1)+selBoxArea(3)))];
                ybounds = [min(selBoxArea(2),(selBoxArea(2)+selBoxArea(4))) max(selBoxArea(2),(selBoxArea(2)+selBoxArea(4)))];
                neurpos = obj.model.neurons_positions;
                selected_objs = 0;
                for i = 1:size(neurpos,1)
                    if neurpos(i,1) > xbounds(1) && neurpos(i,1) < xbounds(2) && neurpos(i,2) > ybounds(1) && neurpos(i,2) < ybounds(2)
%                         selCell = [selCell;obj.model.neuron_objects(i)];
                        selected_objs = selected_objs + 1;
                        viewind = find(~all(sum(neurpos-obj.view.graphic_objects.Neurons(i).Position(1:2)-.5*CanvasConstants.NEURON_size,2),2),1);
                        obj.view.graphic_objects.Neurons(viewind).Selected = 'on';
                    end
                end
                %obj.selection = selCell;
                delete(obj.view.graphic_objects.axes.Children(selInd));
                obj.selection_box = 'off';
                if selected_objs > 0
                    obj.setupContextMenu('groupselected');
                end
%              elseif strcmp(obj.fig.SelectionType,'normal')
            elseif ~isempty(obj.selection) && all(isvalid(obj.selection))
                %if the user left-clicks, clear the selections. Don't clear selections if right-click
                %when jumping from one neuron selection to another, deselect previously selected neuron/s
%                 if size(obj.selection,2) > 1
%                     %obj.clearSelections();
                if ~strcmp(obj.fig.SelectionType,'alt')
                    if size(obj.selection,2) < 2
                        viewNeurPos = cell2mat({obj.view.graphic_objects.Neurons.Position}');
                        selectionPos = obj.selection.Position;
                        [~,deselException] = min(sum(abs(viewNeurPos-selectionPos),2));
                        for i=1:size(viewNeurPos,1)
                            if i~=deselException
                                obj.view.graphic_objects.Neurons(i).Selected = 'off';
                            end
                        end
                    else
                        obj.clearSelections();
                    end
                end
            else
                obj.clearSelections();
            end
            
            obj.fig.WindowButtonMotionFcn = @(~,~) obj.defaultMotion();
            
        end
        %% deleteSelection: deletes item from view and model
        function deleteSelection(obj)
            if ~isempty(obj.selection)
                numSelected = length(obj.selection);
                selTags = {obj.selection.Tag};
                    for ii = 1:numSelected
                        if strcmp(selTags{ii},'n')
                            model_n_pos = obj.model.neurons_positions;
                            sel = obj.selection(ii);
                            pos = sel.Position(1:2) + sel.Position(3:4)/2;
                            [~,model_index] = min(sum(abs(model_n_pos-pos),2));
                        end
                        obj.model.deleteItem(selTags{ii}, model_index);
                    end
                obj.emptySelection();
            end
        end
        %% emptySelection: function for when you click the canvas but don't select an object
        function emptySelection(obj,link_on)
            if nargin == 1
                link_on = 0;
            end
            
            if link_on
                delete(obj.view.graphic_objects.axes.Children(1));
            end  
            
            obj.view.disableForm({'n','l','stimulus'});
            obj.selection = gobjects(0);
            obj.graphic_objects.editMenu.Children(3).Enable = 'off';
            if all(size(obj.view.graphic_objects.Neurons))
                [obj.view.graphic_objects.Neurons.Selected] = deal('off');
            end
        end
        %% setSelection: instantiates the selected object within the controller
        function setSelection(obj, sel)
            obj.selection = sel;
            obj.graphic_objects.editMenu.Children(3).Enable = 'on';
            uistack(obj.selection, 'top');
            obj.selpos = sel.Position;
        end
        %% defaultKeyPress: associates keyboard keys with controller actions
        function defaultKeyPress(obj, data)
            
            if isempty(data.Modifier) && (strcmpi(data.Key, 'delete'))        
                obj.deleteSelection();                
            end
            
        end
        %% modelChanged: sets controller properties to indicate model is unsaved
        function modelChanged(obj)
            
            obj.saved_to_file = false;
            obj.saved_to_workspace = false;
            obj.updateFigureName();
            
        end
        %% updateFigureName: update the figure name to indicate whether model is unsaved
        function updateFigureName(obj)
            
            if obj.saved_to_workspace
                mod1 = '';
            else
                mod1 = '*';
            end
            if obj.has_current_file
                if obj.saved_to_file
                    mod2 = '';
                else
                    mod2 = '*';
                end
                [~,name] = fileparts(obj.current_file);
                cf = [' | ' name mod2];
            else
                cf = '';
            end
            
            obj.fig.Name = ['Canvas Setup' mod1 cf];
            
        end
        %% mOpenModel: remnant of Mars Project, unused
%         function mOpenModel(obj)
%             
%             obj.fig.Pointer = 'watch';
%             drawnow;
%             open_system('SimulationModel')
%             obj.fig.Pointer = 'arrow';
%             
%         end
        %% clearSelections: de-select all active selections
        function clearSelections(obj)
            if ~isempty(obj.view.graphic_objects.Neurons)
                selneurs = sum(contains({obj.view.graphic_objects.Neurons.Selected},'on'));
                if selneurs > 1
                        [obj.view.graphic_objects.Neurons.Selected] = deal('off'); 
                end
            end
        end
    end
end