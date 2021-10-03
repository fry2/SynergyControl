classdef CanvasModel < handle
    
    properties (SetAccess = public, GetAccess=public)
        proj_params = struct()
           
        neurons_positions = zeros(0,2)
        num_neurons = 0
        neuron_objects = struct()
        
        link_ends
        num_links = 0
        link_objects = struct()
        
        muscle_positions = zeros(0,2)
        num_muscles = 0
        muscle_objects = struct()
        
        adapter_positions = zeros(0,2)
        num_adapters = 0
        adapter_objects = struct()
        
        num_stims = 0
        stimulus_objects = struct()
        
        num_datatools = 0;
        datatool_objects = struct()
        
        synapse_types = struct()
        
        obstacles_positions = zeros(0,2)
        num_obstacles = 0        
        
        dtsim = 1e-3;
        tmax = 10;
    end
    
    events (NotifyAccess = private)
        itemAdded
        itemMoved
        itemDeleted
        modelChanged
        modelDeleted
        linkAdded
    end
    
    methods % constructor and destructor
        
        function obj  = CanvasModel()
            obj.setupDefaultSynapseTypes();
            obj.setupProjParams();
        end
        
        function delete(obj)
            obj.notify('modelDeleted');
        end
        
    end
    
    methods (Access = public)
        %% setupSynapseTypes
        function setupDefaultSynapseTypes(obj)
            props = [{'SignalTransmission1','delE',194,'k',1,'max_syn_cond',.115};...
                    {'SignalModulation1','delE',-40,'k',1,'max_syn_cond',.558};...
                    {'SignalModulation2','delE',0,'c',0.05,'max_syn_cond',19};...
                    {'SignalModulation3','delE',-1,'c',0,'max_syn_cond',20}];
            obj.createSynapseType(props);
        end
        %% setupProjParams
        function setupProjParams(obj)
            obj.proj_params.simendtime = 10.01; % times in s
            obj.proj_params.physicstimestep = 0.54; % dt in ms
        end
        %% addItem
        function addItem(obj, type_char, pos, bounds)
            % ADDITEM type_char, pos, bounds
            
            if nargin == 1
                type_char = 'n';
                pos = [50 50];
                bounds = [1000 1000];
            end
            
            %This will change if the randomized position generated in the Controller creates an
            %object that falls outside the bounds of the canvas. If so, it changes the position
            %so that the object doesn't go out of bounds
            
%             pos = CanvasConstants.ConstrainedPosition(type_char, pos);
            %pos = obj.ConstrainedPosition(type_char, pos, bounds);
            
            if type_char == 'n' %neuron
                obj.neurons_positions(end+1,:) = pos;
                obj.num_neurons = size(obj.neurons_positions, 1);
                if obj.num_neurons == 1
                    obj.neuron_objects = obj.create_neuron(pos);
                else
                    obj.neuron_objects(obj.num_neurons,1) = obj.create_neuron(pos);
                end
                index = obj.num_neurons;
            elseif strcmp(type_char,'stimulus')
                [~,modelInd] = min(abs(sum(obj.neurons_positions-pos,2)));
                obj.create_stimulus(modelInd);
                index = modelInd;
            else %obstacle
                obj.obstacles_positions(end+1,:) = pos;
                obj.num_obstacles = size(obj.obstacles_positions, 1);
                index = obj.num_obstacles;
            end
            
            obj.notify('itemAdded', CanvasModelEventData(type_char, index));
            obj.notify('modelChanged');
            
        end
        %% addMuscle
        function addMuscle(obj,muscle_name,linkedID,pos)
            obj.muscle_positions(end+1,:) = pos;
            obj.num_muscles = size(obj.muscle_positions, 1);
            if obj.num_muscles == 1
                obj.muscle_objects = obj.create_muscle(muscle_name,linkedID,pos);
            else
                obj.muscle_objects(obj.num_muscles,1) = obj.create_muscle(muscle_name,linkedID,pos);
            end
            index = obj.num_muscles;
            
            obj.notify('itemAdded', CanvasModelEventData('muscle', index));
            obj.notify('modelChanged');
        end
        %% addAdapter
        function addAdapter(obj,neuron,muscle,pos)
            obj.adapter_positions(end+1,:) = pos;
            obj.num_adapters = size(obj.adapter_positions, 1);
            if obj.num_adapters == 1
                obj.adapter_objects = obj.create_adapter(neuron,muscle,pos);
            else
                obj.adapter_objects(obj.num_adapters,1) = obj.create_adapter(neuron,muscle,pos);
            end
            index = obj.num_adapters;
            
            obj.notify('itemAdded', CanvasModelEventData('adapter', index));
            obj.notify('modelChanged');
        end
        %% addStimulus
        function addStimulus(obj,target,type)
            if isempty(type) || ~any(contains({'tc','dc'},strtrim(lower(type))))
                disp('Stimulus type must be either ''tc'' or ''dc''. No stimulus added.')
                return
            else
                type = strtrim(lower(type));
            end
            try obj.stimulus_objects(1).name;            
                numStims = size(obj.stimulus_objects,1);
            catch
                numStims = 0;
            end 
            
            obj.num_stims = numStims;
            if numStims == 0
                obj.stimulus_objects = obj.create_stimulus(target,type);
            else
                obj.stimulus_objects(numStims+1,1) = obj.create_stimulus(target,type);
            end
            index = numStims;
            
            obj.notify('itemAdded', CanvasModelEventData('stimulus', index));
            obj.notify('modelChanged');
        end
        %% addLink_old
        function addLink_old(obj,start_ind,end_ind,beg,ennd,linktype)
            obj.link_ends(end+1,:) = [beg ennd];
            obj.num_links = size(obj.link_ends, 1);
            index = obj.num_links;
            
            if obj.num_links == 1
                obj.link_objects = obj.create_link([beg ennd]);
            else
                obj.link_objects(obj.num_links,1) = obj.create_link([beg ennd]);
            end
            
            numLinks = obj.num_links;
            
            try obj.link_objects(1).ID;           
                while sum(strcmp({obj.link_objects.ID}, ['link',num2str(numLinks),'-ID']))
                                numLinks = numLinks + 1;
                end
            catch
                numLinks = obj.num_links;
            end
            
            obj.link_objects(obj.num_links,1).ID = ['link',num2str(numLinks),'-ID'];
            obj.link_objects(obj.num_links,1).origin_ID = obj.neuron_objects(start_ind).ID;
            obj.link_objects(obj.num_links,1).destination_ID = obj.neuron_objects(end_ind).ID;
            obj.link_objects(obj.num_links,1).origin_cdata_num = start_ind-1;
            obj.link_objects(obj.num_links,1).destination_cdata_num = end_ind-1;
            obj.link_objects(obj.num_links,1).synaptictype = linktype;
            
            obj.neuron_objects(start_ind).outlinks = [obj.neuron_objects(start_ind).outlinks;obj.link_objects(obj.num_links)];
            obj.neuron_objects(end_ind).inlinks = [obj.neuron_objects(end_ind).inlinks;obj.link_objects(obj.num_links)];
            obj.neuron_objects(start_ind).outlink_IDs = [obj.neuron_objects(start_ind).outlink_IDs;{obj.link_objects(obj.num_links).ID}];
            obj.neuron_objects(end_ind).inlink_IDs = [obj.neuron_objects(end_ind).inlink_IDs;{obj.link_objects(obj.num_links).ID}];
            
            if ~isempty(obj.neuron_objects(start_ind).totmem)
                obj.updateStimModel(end_ind);
            end

            obj.notify('linkAdded', CanvasModelEventData(linktype,[index numLinks],[beg ennd]));
        end
        %% addLink
        function addLink(obj,node1,node2,linktype)
            nodes = {node1, node2};
            for jj = 1:2
                node = nodes{jj};
                switch node.type
                    case 'n'
                        obj_holder = 'neuron_objects';
                    case 'adapter'
                        obj_holder = 'adapter_objects';
                    case 'muscle'
                        obj_holder = 'muscle_objects';
                end
                   link_inds(jj) = find(contains({obj.(obj_holder).name},node.name),1,'first');
                   link_obj_holders{jj} = obj_holder;
            end
            
            if any(ismember(link_obj_holders,'neuron_objects')) && any(ismember(link_obj_holders,'muscle_objects'))
                neuron = nodes{ismember(link_obj_holders,'neuron_objects')};
                muscle = nodes{ismember(link_obj_holders,'muscle_objects')};
                adLoc = [muscle.location(1) (neuron.location(2)+muscle.location(2))/2];
                obj.addAdapter(neuron,muscle,adLoc)
                adInd = size(obj.adapter_objects,1);
                obj.addLink(neuron,obj.adapter_objects(adInd),'adapter')
                obj.addLink(obj.adapter_objects(adInd),muscle,'adapter')
                return
            end
            
            obj.link_ends(end+1,:) = [node1.location node2.location];
            obj.num_links = size(obj.link_ends, 1);
            numLinks = obj.num_links;
            
            if numLinks == 1
                obj.link_objects = obj.create_link([node1.location node2.location]);
            else
                obj.link_objects(numLinks,1) = obj.create_link([node1.location node2.location]);
            end
            
            try obj.link_objects(1).ID;           
                while sum(strcmp({obj.link_objects.ID}, ['link',num2str(numLinks),'-ID']))
                                numLinks = numLinks + 1;
                end
            catch
                numLinks = obj.num_links;
            end
            
            obj.link_objects(numLinks,1).ID = ['link',num2str(numLinks),'-ID'];
            obj.link_objects(numLinks,1).origin_ID = node1.ID;
            obj.link_objects(numLinks,1).destination_ID = node2.ID;
            %%obj.link_objects(numLinks,1).origin_cdata_num = start_ind-1;
            %%obj.link_objects(numLinks,1).destination_cdata_num = end_ind-1;
            obj.link_objects(numLinks,1).synaptictype = linktype;
            
            obj.(link_obj_holders{1})(link_inds(1)).outlinks = [obj.(link_obj_holders{1})(link_inds(1)).outlinks;obj.link_objects(obj.num_links)];
            obj.(link_obj_holders{2})(link_inds(2)).inlinks = [obj.(link_obj_holders{2})(link_inds(2)).inlinks;obj.link_objects(obj.num_links)];
            obj.(link_obj_holders{1})(link_inds(1)).outlink_IDs = [obj.(link_obj_holders{1})(link_inds(1)).outlink_IDs;{obj.link_objects(obj.num_links).ID}];
            obj.(link_obj_holders{2})(link_inds(2)).inlink_IDs = [obj.(link_obj_holders{2})(link_inds(2)).inlink_IDs;{obj.link_objects(obj.num_links).ID}];
            
%             obj.neuron_objects(start_ind).outlinks = [obj.neuron_objects(start_ind).outlinks;obj.link_objects(obj.num_links)];
%             obj.neuron_objects(end_ind).inlinks = [obj.neuron_objects(end_ind).inlinks;obj.link_objects(obj.num_links)];
%             obj.neuron_objects(start_ind).outlink_IDs = [obj.neuron_objects(start_ind).outlink_IDs;{obj.link_objects(obj.num_links).ID}];
%             obj.neuron_objects(end_ind).inlink_IDs = [obj.neuron_objects(end_ind).inlink_IDs;{obj.link_objects(obj.num_links).ID}];
            
%             if ~isempty(obj.neuron_objects(start_ind).totmem)
%                 obj.updateStimModel(end_ind);
%             end

            obj.notify('linkAdded', CanvasModelEventData(linktype,[numLinks numLinks],[node1.location node2.location]));
        end
        %% addDatatool
        function addDatatool(obj,inParams)
            parCell = iscell(inParams);
            if parCell
                name = inParams{1};
                while ~ischar(name)
                    name = input('First input parameter must be a string name for the datatool.\nPlase enter a string for the datatool name now: \n');
                end
            else
                name = char(inParams);
            end
            
            if contains(name,' ')
                warning('Datatool name contains a space, removing spaces.')
                name = name(~isspace(name));
            end
            
            try obj.datatool_objects(1).name;            
                numDatatools = size(obj.datatool_objects,1);
            catch
                numDatatools = 0;
            end 
            
            obj.num_datatools = numDatatools;
            
            if numDatatools == 0
                datInd = 1;
                obj.datatool_objects = obj.create_datatool(name);
            else
                if any(contains({obj.datatool_objects.name},name))
                    warning(['CanvasModel.addDatatool: Datatool with name ''',name,''' already exists. Overwriting datatool with new information.'])
                    datInd = find(contains({obj.datatool_objects.name},name),1,'first');
                    obj.datatool_objects(datInd,1) = obj.create_datatool(name);
                else
                    datInd = numDatatools+1;
                    obj.datatool_objects(datInd,1) = obj.create_datatool(name);
                end
            end
            
            if parCell
                for ii = 2:2:length(inParams)
                    parName = char(inParams{ii});
                    parVal = double(inParams{ii+1});
                    if any(contains(fields(obj.datatool_objects),parName))
                        obj.datatool_objects(datInd,1).(parName) = parVal;
                    else
                        warning(['CanvasModel.addDatatool: ',parName,' for datatool ',name,' does not exist. Information not added to CanvasModel.'])
                    end
                end
            end
            
            obj.notify('itemAdded', CanvasModelEventData('datatool', numDatatools));
            obj.notify('modelChanged');
        end
        %% addDTaxes
        function addDTaxes(obj,datatool,target,datatype)
            if ~any(contains({'MembraneVoltage','Tension','Activation','Tl','MuscleLength'},datatype))
                warning([target.name,' datatype (',datatype,') incorrect in addDTaxes'])
            end
            
            if ischar(datatool) || isstring(datatool)
                % Datatool has been specified by a name string
                dind = find(contains({obj.datatool_objects.name},datatool),1,'first');
                if isempty(dind)
                    error(['Specified datatool ',datatool,' does not exist.'])
                else
                    datatool = obj.datatool_objects(dind);
                end
            end
            
            try datatool.axes_objects(1).name;            
                numDTaxes = size(datatool.axes_objects,1);
            catch
                numDTaxes = 0;
            end
            
            if numDTaxes == 0
                datatool.axes_objects = obj.create_dtaxes(target,datatype);
            else
                datatool.axes_objects(numDTaxes+1,1) = obj.create_dtaxes(target,datatype);
            end
            
            dtInd = find(ismember({obj.datatool_objects.name},datatool.name));
            obj.datatool_objects(dtInd,1) = datatool;
            
            obj.notify('itemAdded', CanvasModelEventData('dtaxes', numDTaxes));
            obj.notify('modelChanged');
        end
        %% moveItem
        function moveItem(obj, type_char, index, pos, bounds)
            
            if size(index,1) > 1
                index = index(1);
            end
            
            modelneurons = {obj.neuron_objects(:).ID};
            
%             pos = CanvasConstants.ConstrainedPosition(type_char, pos);
            pos = obj.ConstrainedPosition(type_char, pos, bounds);
            
            if type_char == 'n' %neuron
                obj.neurons_positions(index, :) = pos;
                obj.neuron_objects(index).location = pos;
                if ~isempty(obj.neuron_objects(index).inlinks) || ~isempty(obj.neuron_objects(index).outlinks)
                    modellinks = {obj.link_objects(:).ID};
                    if isempty(obj.neuron_objects(index).inlinks)
                        numInlinks = 0;
                    else
                        numInlinks = size(obj.neuron_objects(index).inlinks,1);
                        inlinkcell = obj.neuron_objects(index).inlink_IDs;
                    end
                    
                    if isempty(obj.neuron_objects(index).outlinks)
                        numOutlinks = 0;
                    else
                        numOutlinks = size(obj.neuron_objects(index).outlinks,1);
                        outlinkcell = obj.neuron_objects(index).outlink_IDs;
                    end
                            
                    for i = 1:numInlinks
                        % Finding where everything is (indices are different for different lists)
                        linkID = obj.neuron_objects(index).inlink_IDs{i};
                        modellinkInd = find(contains(modellinks,linkID),1);
                        inlinkInd = find(strcmp(inlinkcell, linkID)==1);
                        proxindex = find(contains(modelneurons,obj.neuron_objects(index).inlinks(inlinkInd).origin_ID),1);
                        proxoutlinkcell = {obj.neuron_objects(proxindex).outlinks.ID};
                        proxoutlinkInd = strcmp(proxoutlinkcell, linkID)==1;
                        
                        % Assigning revised position to each object
                        obj.link_objects(modellinkInd).end = pos;
                        obj.link_ends(modellinkInd,3:4) = pos;
                        obj.neuron_objects(index).inlinks(inlinkInd).end = pos;
                        obj.neuron_objects(proxindex).outlinks(proxoutlinkInd).end = pos;
                    end
                    
                    for j = 1:numOutlinks
                        % Finding where everything is (indices are different for different lists)
                        linkID = obj.neuron_objects(index).outlink_IDs{j};
                        modellinkInd = find(contains(modellinks,linkID),1);
                        outlinkInd = find(strcmp(outlinkcell, linkID)==1);
                        distalindex = find(contains(modelneurons,obj.neuron_objects(index).outlinks(outlinkInd).destination_ID),1);
                        distalinlinkcell = {obj.neuron_objects(distalindex).inlinks.ID};
                        distalinlinkInd = strcmp(distalinlinkcell, linkID)==1;
                        
                        % Assigning revised position to each object
                        obj.link_objects(modellinkInd).start = pos;
                        obj.link_ends(modellinkInd,1:2) = pos;
                        obj.neuron_objects(index).outlinks(outlinkInd).start = pos;
                        obj.neuron_objects(distalindex).inlinks(distalinlinkInd).start = pos;
                    end
                end
            else %obstacle
                obj.obstacles_positions(index, :) = pos;
            end
            
            obj.notify('itemMoved', CanvasModelEventData(type_char,index));
            obj.notify('modelChanged');
            
        end
        %% deleteItem
        function deleteItem(obj, type_char, index)
            
            obj.notify('itemDeleted', CanvasModelEventData(type_char, index));
            
            if type_char == 'n' %neuron
                if ~isempty(obj.neuron_objects(index).outlink_IDs) || ~isempty(obj.neuron_objects(index).inlink_IDs)
                    if size(obj.neuron_objects(index).inlinks,2) == 0
                        numInlinks = 0;
                    else
                        numInlinks = size(obj.neuron_objects(index).inlink_IDs,1);
                    end
                    
                    if size(obj.neuron_objects(index).outlinks,2) == 0
                        numOutlinks = 0;
                    else
                        numOutlinks = size(obj.neuron_objects(index).outlink_IDs,1);
                    end
                    while numOutlinks > 0
                        linkIDs = {obj.link_objects.ID};
                        linkIDstring = obj.neuron_objects(index).outlinks(1).ID;
                        linkind = find(strcmp(linkIDs, linkIDstring)==1);
                        dNeurInd = find(contains({obj.neuron_objects.ID},obj.link_objects(linkind).destination_ID));
                        obj.deleteItem('l', linkind);
                        obj.updateStimModel(dNeurInd);
                        numOutlinks = numOutlinks - 1;
                    end
                    while numInlinks > 0
                        linkIDs = {obj.link_objects.ID};
                        linkIDstring = obj.neuron_objects(index).inlinks(1).ID;
                        linkind = find(strcmp(linkIDs, linkIDstring)==1);
                        obj.deleteItem('l', linkind);
                        numInlinks = numInlinks - 1;
                    end
                end
                if size(obj.neurons_positions,1) == 1
                    obj.neuron_objects = [];
                    obj.neurons_positions = [];
                else
                    obj.neuron_objects(index) = [];
                    obj.neurons_positions(index, :) = [];
                end
                obj.num_neurons = size(obj.neurons_positions, 1);
                obj.update_link_cdata();
            elseif type_char == 'l' %link
                link = obj.link_objects(index);
                neuronIDs = {obj.neuron_objects.ID};
                
                startind = find(strcmp(neuronIDs, link.origin_ID)==1);
                outind = find(strcmp({obj.neuron_objects(startind).outlinks.ID}, link.ID)==1);
                
                if size(obj.neuron_objects(startind).outlinks,1) == 1
                    obj.neuron_objects(startind).outlinks = [];
                    obj.neuron_objects(startind).outlink_IDs = [];
                else
                    obj.neuron_objects(startind).outlinks(outind) = [];
                    obj.neuron_objects(startind).outlink_IDs(outind) = [];
                end

                endind = find(strcmp(neuronIDs, link.destination_ID)==1);
                inind = find(strcmp({obj.neuron_objects(endind).inlinks.ID}, link.ID)==1);
                if size(obj.neuron_objects(endind).inlinks,1) == 1
                    obj.neuron_objects(endind).inlinks = [];
                    obj.neuron_objects(endind).inlink_IDs = [];
                else
                    obj.neuron_objects(endind).inlinks(inind) = [];
                    obj.neuron_objects(endind).inlink_IDs(inind) = [];
                end
                    
                if size(obj.link_objects,1) == 1
                    obj.link_objects = [];
                    obj.link_ends = [];
                else
                    obj.link_objects(index) = [];
                    obj.link_ends(index,:) = [];
                end
                obj.num_links = size(obj.link_ends,1);
            else %obstacle
                obj.obstacles_positions(index, :) = [];
                obj.num_obstacles = size(obj.obstacles_positions, 1);
            end
            
            
            obj.notify('modelChanged');
            
        end
        %% getData
        function S = getData(obj)

            for i=1:size(obj.link_objects,1)
                S.Link_Objects(i) = obj.link_objects(i);
            end
            
            for i=1:size(obj.neuron_objects,1)
                S.Neuron_Objects(i) = obj.neuron_objects(i);
            end
            
        end
        %% setData
        function setData(obj, S)
            
            for k=obj.num_neurons:-1:1
                obj.deleteItem('n', k);
            end
            
            for k=obj.num_obstacles:-1:1
                obj.deleteItem('o', k);
            end
            
            for k=1:size(S.NeuronsPositions, 1)
                obj.addItem('n', S.NeuronsPositions(k,:));
            end
            
            for k=1:size(S.ObstaclesPositions, 1)
                obj.addItem('o', S.ObstaclesPositions(k,1:2) + ...
                    S.ObstaclesPositions(k,3:4)/2);
            end
            
        end
        %% newCanvas
        function newCanvas(obj, numAddNeurons ,bounds)
            
            for k=obj.num_neurons:-1:1
                obj.deleteItem('n', k);
            end
            
            for k=1:numAddNeurons
                obj.addItem('n', [CanvasConstants.NEURON_size(1) CanvasConstants.NEURON_size(2)*1.5*k],bounds);
            end
            
        end
        %% create_link
        function link = create_link(~,pos)
            link = struct;
            link.ID = '';
            link.origin_cdata_num = [];
            link.destination_cdata_num = [];
            link.start = pos(1:2);
            link.end = pos(3:4);
            link.origin_ID = '';
            link.origin_cdata_num = [];
            link.destination_ID = '';
            link.destination_cdata_num = [];
            link.synaptictype = {};
            link.assemblyfile = 'IntegrateFireGUI.dll';
            link.behavior = 'IntegrateFireGUI.DataObjects.Behavior.Synapse';
        end
        %% create_synapseType
        function createSynapseType(obj,props)
            transnames = {'transmission','trans','hyperpolarizing'};
            modnames = {'modulation','mod','depolarizing'};
            if isempty(fields(obj.synapse_types))
                numSynapseTypes = 0;
            else
                numSynapseTypes = size(obj.synapse_types,1);
            end
            for j=1:size(props,1)
                synaptic_properties = props(j,:);
                synNum = numSynapseTypes+j;
                obj.synapse_types(synNum,1).ID = [lower(synaptic_properties{1,1}),'-ID'];
                for k = 2:2:size(props,2)
                    obj.synapse_types(synNum,1).name = props{j,1};
                    property_string = props{j,k};
                    property_value = props{j,k+1};
                    obj.synapse_types(synNum,1).(property_string) = property_value;
                end
                    obj.synapse_types(synNum,1).presyn_thresh = -60;
                    obj.synapse_types(synNum,1).presyn_sat = -40;
                    param_ind = logical([0,strcmp(synaptic_properties(1:end-1),'delE')]);
                    delE = synaptic_properties{param_ind};
                    obj.synapse_types(synNum,1).equil_pot = delE-60;
                if ~any(strcmp(synaptic_properties,'max_syn_cond'))
                    if any(strcmp(synaptic_properties,'delE'))
                        if any(strcmp(synaptic_properties,'k'))
                            param_ind = logical([0,strcmp(synaptic_properties(1:end-1),'k')]);
                            mod_param = synaptic_properties{param_ind};
                        else any(strcmp(synaptic_properties,'c'))
                            param_ind = logical([0,strcmp(synaptic_properties(1:end-1),'c')]);
                            mod_param = synaptic_properties{param_ind};
                        end
                    else
                        delE = input('No delE was provided. Please enter a delE value in mV:\n');
                    end
                    tempSynCond = (mod_param*20)/(delE-mod_param*20);
                    if log10(tempSynCond) < -3
                        % If the calculated maximum synaptic conductance is very small, just set it to 0
                        % This is done primarily to 1) avoid Animatlab errors and 2) deal with constant joint angle inputs
                        tempSynCond = 0;
                    end
                    obj.synapse_types(synNum,1).max_syn_cond = tempSynCond;
                end
                if contains(lower(synaptic_properties{1,1}),modnames)
                    obj.synapse_types(synNum,1).arrow_dest_style = 'Circle';
                elseif contains(lower(synaptic_properties{1,1}),transnames)
                    obj.synapse_types(synNum,1).arrow_dest_style = 'Fork';
                else
                    obj.synapse_types(synNum,1).arrow_dest_style = 'Fork';
                end
                    obj.synapse_types(synNum,1).arrow_dest_size = 'Small';
                    obj.synapse_types(synNum,1).arrow_dest_angle = 'deg30';
                    obj.synapse_types(synNum,1).arrow_dest_filled = 'False';
                    obj.synapse_types(synNum,1).arrow_mid_style = 'None';
                    obj.synapse_types(synNum,1).arrow_mid_size = 'Small';
                    obj.synapse_types(synNum,1).arrow_mid_angle = 'deg30';
                    obj.synapse_types(synNum,1).arrow_mid_filled = 'False';
                    obj.synapse_types(synNum,1).arrow_origin_style = 'None';
                    obj.synapse_types(synNum,1).arrow_origin_size = 'Small';
                    obj.synapse_types(synNum,1).arrow_origin_angle = 'deg30';
                    obj.synapse_types(synNum,1).arrow_origin_filled = 'False';
            end
        end
        %% create_neuron
        function neuron = create_neuron(obj,pos)
            neuron = struct;
            numNeurons = obj.num_neurons;
            try obj.neuron_objects(1).name;            
                while sum(strcmp({obj.neuron_objects.name}, ['neur',num2str(numNeurons)]))
                                numNeurons = numNeurons + 1;
                end
            catch
                numNeurons = obj.num_neurons;
            end
            neuron.name = strcat('neur',num2str(numNeurons));
            neuron.ID = [neuron.name,'-ID'];
            neuron.location = pos;
            neuron.nsize = CanvasConstants.NEURON_size;
            neuron.outlinks = [];
            neuron.outlink_IDs = {};
            neuron.inlinks = [];
            neuron.inlink_IDs = {};
            neuron.enabled = CanvasConstants.NEURON_enabled;
            neuron.restingpotential = CanvasConstants.NEURON_restingpotential;
            neuron.timeconstant = CanvasConstants.NEURON_timeconstant;
            neuron.initialthreshold = CanvasConstants.NEURON_initialthreshold;
            neuron.relativeaccomodation = CanvasConstants.NEURON_relativeaccomodation;
            neuron.accomodationtimeconstant = CanvasConstants.NEURON_accomodationtimeconstant;
            neuron.AHPconductance = CanvasConstants.NEURON_AHPconductance;
            neuron.AHPtimeconstant = CanvasConstants.NEURON_AHPtimeconstant;
            neuron.ca_act_ID = [neuron.ID,'-act'];
            neuron.ca_act_midpoint = CanvasConstants.NEURON_ca_act_midpoint;
            neuron.ca_act_slope = CanvasConstants.NEURON_ca_act_slope;
            neuron.ca_act_timeconstant = CanvasConstants.NEURON_ca_act_timeconstant;
            neuron.ca_deact_ID = [neuron.ID,'-deact'];
            neuron.ca_deact_midpoint = CanvasConstants.NEURON_ca_deact_midpoint;
            neuron.ca_deact_slope = CanvasConstants.NEURON_ca_deact_slope;
            neuron.ca_deact_timeconstant =CanvasConstants.NEURON_ca_deact_timeconstant;
            neuron.tonicstimulus = CanvasConstants.NEURON_tonicstimulus;
            neuron.tonicnoise = CanvasConstants.NEURON_tonicnoise;
            neuron.type = 'n';
            neuron.stimulus = {};
            neuron.totmem = {};
            neuron.color = -7876870;
        end
        %% create_muscle
        function muscle = create_muscle(obj,muscle_name,linkedID,pos)
            muscle = struct;
            numMuscles = obj.num_muscles;
            try obj.muscle_objects(1).name;            
                while sum(strcmp({obj.muscle_objects.name}, ['musc',num2str(numMuscles)]))
                                numMuscles = numMuscles + 1;
                end
            catch
                numMuscles = obj.num_muscles;
            end
            muscle.name = muscle_name;
            muscle.ID = [muscle_name,'-ID'];
            muscle.location = pos;
            muscle.linkedID = linkedID;
            muscle.inlinks = [];
            muscle.inlink_IDs = {};
            muscle.outlink_IDs = {};
            muscle.outlinks = [];
            muscle.type = 'muscle';
            muscle.size = CanvasConstants.MUSCLE_size;
        end
        %% create_adapter
        function adapter = create_adapter(obj,neuron,muscle,pos)
            adapter = struct;
            numAdapters = obj.num_adapters;
            try obj.adapter_objects(1).name;            
                while sum(strcmp({obj.adapter_objects.name}, ['musc',num2str(numAdapters)]))
                                numAdapters = numAdapters + 1;
                end
            catch
                numAdapters = obj.num_adapters;
            end          
            adapter.name = ['ad-',muscle.name(1:end-7)];
            adapter.ID = [adapter.name,'-ID'];
            adapter.location = pos;
            adapter.size = CanvasConstants.ADAPTER_size;
            adapter.origin_node_ID = neuron.ID;
            adapter.origin_node_name = neuron.name;
            adapter.destination_node_ID = muscle.ID;
            adapter.destination_linked_ID = muscle.linkedID;
            adapter.origin_node_name = muscle.name;
            adapter.inlinks = [];
            adapter.outlinks = [];
            adapter.inlink_IDs = {};
            adapter.outlink_IDs = [];
            adapter.type = CanvasConstants.ADAPTER_type;
            adapter.gain_profile_ID = [adapter.ID,'-',CanvasConstants.ADAPTER_gain_profile_ID];
            
            neurInd = find(contains({obj.neuron_objects.name},neuron.name));
            muscInd = find(contains({obj.muscle_objects.name},muscle.name));
            %obj.addLink(neurInd,muscInd,neuron.location,muscle.location,linktype)
        end
        %% create_stimulus
        function stimulus = create_stimulus(obj,target,type)
            if strcmp(target.type,'n')
                stimulus = struct; 
                % To allow multiple stimuli to be added to a neuron, we need to 1) check that the proposed stim is a duplicate,
                % 2) if so, increment a counter until the new stimulus is unique
                stimDup = 1;
                try obj.stimulus_objects(1).name;            
                    temp_name = ['st',upper(type),num2str(stimDup),'-',target.name];
                    while any(contains({obj.stimulus_objects.name},temp_name))
                        stimDup = stimDup + 1;
                        temp_name = ['st',upper(type),num2str(stimDup),'-',target.name];
                    end
                catch
                    temp_name = ['st',upper(type),num2str(stimDup),'-',target.name];
                end 

                stimulus.name = temp_name;
                stimulus.ID = [stimulus.name,'-ID'];
                stimulus.target_ID = target.ID;
                stimulus.enabled ='True';
                stimulus.starttime = 0;
                stimulus.endtime = obj.proj_params.simendtime;
                stimulus.magnitude = 10;
                stimulus.simeq = [];
                stimulus.projeq = [];
                stimulus.rest_potential = -60;
                stimulus.conductance = 100;
                stimulus.current_wave = [];
                stimulus.current_data_file = [];
                stimulus.muscle_ID = [];
                stimulus.type = strtrim(lower(type));
                
                if strcmp(type,'dc')
                    % In order for a direct current stimulus to work, it needs a muscle ID
                    % Because we have bypassed other internal processes, this can be *any* muscle ID
                    % For simplicity, we just add the first muscle ID in the system.
                    try 
                        obj.muscle_objects(1).name;
                        stimulus.muscle_ID = obj.muscle_objects(1).linkedID;
                    catch ME
                        error(['Trying to add a direct current stimulus without any muscles in the object. You must have at least one muscle in the object'...
                            ' in order for the direct current stimulus to work properly.'])
                    end
                else
                    stimulus.muscle_ID = [];
                end
            else
                disp('Stims only for neurons. Please select a neuron.')
            end
        end
        %% create_datatool
        function datatool = create_datatool(~,name)
           name = char(name);
           datatool = struct;
           datatool.name = name;
           datatool.ID = ['dt-',name,'-ID'];
           datatool.tfID = ['dt-',name,'-tfID'];
           datatool.starttime = 0;
           datatool.endtime = 10;
           datatool.collectdatainterval = 0.54; % collect data interval in ms
           datatool.axes_objects = struct();
        end
        %% create_dtaxes
        function dtaxes = create_dtaxes(~,target,datatype)
            dtaxes = struct;
            dtaxes.name = ['dtax-',target.name];
            dtaxes.tdata_ID = [target.name,'-tdata-ID'];
            dtaxes.target_name = target.name;
            if strcmp(target.type,'muscle')
                % Dark red
                dtaxes.linecolor = -6139310;
                dtaxes.target_type = 'muscle';
                dtaxes.target_ID = target.linkedID;
                if isempty(datatype) || ~ischar(datatype)
                    dtaxes.datatype = 'Tension';
                else
                    dtaxes.datatype = datatype;
                end
            else
                dtaxes.target_type = 'neuron';
                dtaxes.linecolor = target.color;
                dtaxes.target_ID = target.ID;
                if isempty(datatype) || ~ischar(datatype)
                    dtaxes.datatype = 'MembraneVoltage';
                else
                    dtaxes.datatype = datatype;
                end
            end
        end
        %% create_animatlab project
        function create_animatlab_project(obj,proj_file)
            revised_file = strcat(proj_file(1:end-6),'_fake.aproj');

            if size(revised_file,2) <= 7
                return
            end
            original_text = importdata(proj_file);
            modified_text = original_text;
            
            %%%Overwrite Project Parameters
            for ii = 1:size(fields(obj.proj_params),1)
                pFields = fields(obj.proj_params);
                pInd = find(contains(lower(original_text),pFields{ii}));
                quoteLocs = strfind(original_text{pInd},'"');
                if ~isempty(quoteLocs)
                    sTemp = original_text{pInd};
                    qTemp = reshape(quoteLocs,[2 3])';
                    scaleStr = sTemp(qTemp(2,1)+1:qTemp(2,2)-1);
                    switch scaleStr
                        case 'None'
                            scaler = 1;
                        case 'milli'
                            scaler = .001;
                        case 'nano'
                            scaler = 1e-9;
                        case 'micro'
                            scaler = 1e-6;
                    end
                    if strcmp(pFields{ii},'physicstimestep')
                        % Timestep will be given in milliseconds always, overwrite the original
                        scaleStr = 'milli';
                        scaler = 1e-3;
                    end
                    pVal = obj.proj_params.(pFields{ii});
                    modified_text{pInd} = [sTemp(1:qTemp(1,1)),num2str(pVal),sTemp(qTemp(1,2):qTemp(2,1)),scaleStr,...
                                            sTemp(qTemp(2,2):qTemp(3,1)),num2str(scaler*pVal),sTemp(qTemp(3,2):end)];
                else
                    sTemp = extractBetween(string(original_text{pInd}),'>','<');
                    if isnan(str2double(sTemp))
                        keyboard
                    else
                        keyboard
                    end
                end
            end
            
            %%%Overwrite the <NervousSystem> neural subsystem information
            if isempty(find(contains(original_text,'<NeuralModules>'),1))
                nervoussystem_inject_start = find(contains(original_text,'<NeuralModules/>'))-1;
            else
                nervoussystem_inject_start = find(contains(original_text,'<NeuralModules>'))-1;
            end
            nervoussystem_inject_end = find(contains(original_text,'</NervousSystem>'));
            [ns_nervoussystem_text,ns_tab_text] = CanvasText(obj).build_ns_text('project');            
            modified_text = [modified_text(1:nervoussystem_inject_start,1);...
                ns_nervoussystem_text;...
                modified_text(nervoussystem_inject_end:end,1)];
            
            %%%Inject Datatool Code
            try obj.datatool_objects(1).name;
                numDatatools = size(obj.datatool_objects,1);
                aform_dt = obj.datatool_objects(1).collectdatainterval;
            catch
                numDatatools = 0;
                aform_dt = obj.proj_params.physicstimestep;
            end
            
            if numDatatools > 0 
                if isempty(find(contains(modified_text,'<ToolViewers/>'),1))
                    datatool_inject = find(contains(modified_text,'</ToolViewers>'));
                    datatool_text = {};
                else
                    datatool_inject = find(contains(modified_text,'<ToolViewers/>'));
                    datatool_text = {'<ToolViewers>'};
                end
                
                for i = 1:numDatatools
                    datatool = obj.datatool_objects(i);
                    if isempty(find(contains(modified_text,['<Name>',datatool.name,'</Name>']),1))
                        datatool_holder = CanvasText(obj).build_datatool(datatool);
                        datatool_text = [datatool_text;datatool_holder];
                        aform_path = [fileparts(proj_file),'\',datatool.name,'.aform'];
                        aform_text = CanvasText(obj).build_aform_text(datatool,'project');
                            fileID = fopen(aform_path,'w');
                            fprintf(fileID,'%s\n',aform_text{:});
                            fclose(fileID);
                    end
                end
                datatool_text{end+1,1} = '</ToolViewers>';
                modified_text = [modified_text(1:datatool_inject-1);...
                                datatool_text;...
                                modified_text(datatool_inject+1:end)];
            end
            
            %%%Modify existing datatools such that they terminate before the end of the simulation
            
            extDTInds = [find(contains(original_text,'<Name>JointMotion')),find(contains(original_text,'<Name>PassiveTension'))];
            extDTs = cell(1,length(extDTInds));
            for ii = 1:length(extDTInds)
                extDTs{ii} = char(extractBetween(string(original_text{extDTInds(ii)}),'>','</'));
            end
            
            %extDTs = {'JointMotion';'PassiveTension'};
            parCell = {'<CollectEndTime',obj.proj_params.simendtime-.01;...
                       '<CollectDataInterval',aform_dt/1000};
            for ii = 1:length(extDTs)
                aform_path = [fileparts(proj_file),'\',extDTs{ii},'.aform'];
                aform_text = importdata(aform_path);
                for jj = 1:size(parCell,1)
                    pInd = find(contains(aform_text,parCell{jj,1}),1,'first');
                    sTemp = aform_text{pInd};
                    qTemp = reshape(strfind(sTemp,'"'),[2 3])';
                    scaler = 1;
                    pVal = parCell{jj,2};
                    aform_text{pInd} = [sTemp(1:qTemp(1,1)),num2str(pVal),sTemp(qTemp(1,2):qTemp(2,1)),'None',...
                                            sTemp(qTemp(2,2):qTemp(3,1)),num2str(scaler*pVal),sTemp(qTemp(3,2):end)];
                end
                fileID = fopen(aform_path,'w');
                fprintf(fileID,'%s\n',aform_text{:});
                fclose(fileID);
                clear fileID aform_text stemp qTemp pVal aform_path
            end
            
            %%%Inject Stimulus Code
            try obj.stimulus_objects(1).name;            
                numStims = size(obj.stimulus_objects,1);
            catch
                numStims = 0;
            end 
            
            % Disable any existing stimuli
            extStims = find(contains(modified_text,'<Stimulus>'));
            for ii = 1:length(extStims)
                enInd = find(contains(modified_text(extStims(ii):end),'<Enabled>'),1,'first')+extStims(ii)-1;
                modified_text{enInd} = '<Enabled>False</Enabled>';
            end
            
            if numStims > 0
                stimuli_inject = find(contains(modified_text,'</Stimuli>'));    
                stim_text = {};
                for i=1:numStims
                    stim_holder = CanvasText(obj).build_stimulus(obj.stimulus_objects(i),'project');
                    stim_text = [stim_text;stim_holder];
                end
                modified_text = [modified_text(1:stimuli_inject-1);...
                                stim_text;...
                                modified_text(stimuli_inject:end)];
            end
            
            %%%Overwrite the <TabbedGroupsConfig> Information
            % First, find the existing DTs that you want to preserve and append them to the ns_tab_text
            for ii = 1:length(extDTs)
                dt1 = find(contains(modified_text,['Page Title="',extDTs{ii},'"']),1);
                dt2 = find(contains(modified_text,['Page Title="',extDTs{ii},'"']),1)+12;
                snip = modified_text(dt1:dt2);
                ns_tab_text = [ns_tab_text;snip];
            end
            
            tab_inject_start = find(contains(modified_text,'&lt;SubSystemID&gt'),1,'first')-9;
            tab_inject_end = find(contains(modified_text,'&lt;/Leaf&gt;'),1,'first');
            modified_text = [modified_text(1:tab_inject_start-1,1);...
                            ns_tab_text;...
                            modified_text(tab_inject_end:end,1)];
            leaf_num = size(find(contains(modified_text,'&lt;Page Title="')),1);
            leafInd = find(contains(modified_text,'&lt;Leaf Count'),1,'first');
            modified_text{leafInd} = ['&lt;Leaf Count="',num2str(leaf_num),'" Unique="7" Space="100"&gt;'];
            
            [~,projName,projExt] = fileparts(revised_file);
            disp(['Project file ',projName,projExt,' has been updated.'])
            fileID = fopen(revised_file,'w');
            fprintf(fileID,'%s\n',modified_text{:});
            fclose(fileID);
        end
        %% create_animatlab simulation
        function create_animatlab_simulation(obj,sim_file)

            original_text = importdata(sim_file);
            modified_text = original_text;
            
            %%%Overwrite Project Parameters
            pFields = fields(obj.proj_params);
            for ii = 1:size(pFields,1)
                pInd = find(contains(lower(original_text),pFields{ii}));
                sTemp = string(original_text{pInd});
                sTemp = char(eraseBetween(sTemp,'>','<'));
                gt = strfind(sTemp,'>');
                if strcmp(pFields{ii},'physicstimestep')
                    pStr = num2str(obj.proj_params.(pFields{ii})/1000);
                else
                    pStr = num2str(obj.proj_params.(pFields{ii}));
                end
                modified_text{pInd} = [sTemp(1:gt(1)),pStr,sTemp(gt(1)+1:end)];
            end
            
            %%%Overwrite the <NervousSystem> neural subsystem information
            if isempty(find(contains(original_text,'<NeuralModules>'),1))
                nervoussystem_inject_start = find(contains(original_text,'<NeuralModules/>'))-1;
            else
                nervoussystem_inject_start = find(contains(original_text,'<NeuralModules>'))-1;
            end
            nervoussystem_inject_end = find(contains(original_text,'</NervousSystem>'));
            ns_nervoussystem_text = CanvasText(obj).build_ns_text('simulation');            
            modified_text = [modified_text(1:nervoussystem_inject_start,1);...
                ns_nervoussystem_text;...
                modified_text(nervoussystem_inject_end:end,1)];
            
            %%%Inject Datatool Code
            try obj.datatool_objects(1).name;
                numDatatools = size(obj.datatool_objects,1);
                aform_dt = obj.datatool_objects(1).collectdatainterval;
            catch
                numDatatools = 0;
                aform_dt = obj.proj_params.physicstimestep;
            end
            
            %%%Modify existing datatools such that they terminate before the end of the simulation
            
            extDTs = {'<Name>JointMotion';...
                      '<Name>PassiveTension'};
            parCell = {'<EndTime',obj.proj_params.simendtime-.01;...
                       '<CollectInterval',aform_dt/1000};
            temp_extDatHolder = {};
            for ii = 1:length(extDTs)
                for jj = 1:size(parCell,1)
                    dtInd = find(contains(modified_text,extDTs{ii}));
                    pInd = find(contains(modified_text(dtInd:end),parCell{jj,1}),1,'first')+dtInd-1;
                    sTemp = string(modified_text{pInd});
                    sTemp = char(eraseBetween(sTemp,'>','<'));
                    gt = strfind(sTemp,'>');
                    pStr = num2str(parCell{jj,2});
                    modified_text{pInd} = [sTemp(1:gt(1)),pStr,sTemp(gt(1)+1:end)];
                end
                dChartEnd = find(contains(modified_text(dtInd-4:end),'</DataChart>'),1,'first')+dtInd-5;
                temp_extDatHolder = [temp_extDatHolder;modified_text((dtInd-4):dChartEnd)];
            end
            
            if ~isempty(extDTs)
                modified_text = [modified_text(1:find(contains(modified_text,'<DataCharts>')));...
                                temp_extDatHolder;...
                                modified_text(find(contains(modified_text,'</DataCharts>')):end)];
            end
            
            
            if numDatatools > 0
                datatool_inject = find(contains(modified_text,'<DataCharts/>'),1);
                if isempty(datatool_inject)
                    datatool_inject0 = find(contains(modified_text,'</DataChart>'),1,'last');
                    datatool_inject1 = find(contains(modified_text,'</DataCharts>'));
                    datatool_text = {};
                else
                    datatool_inject0 = datatool_inject-1;
                    datatool_inject1 = datatool_inject+1;
                    datatool_text = {'<DataCharts>'};
                end
                
                for i = 1:numDatatools
                    datatool = obj.datatool_objects(i);
                    chart_holder = CanvasText(obj).build_aform_text(datatool,'simulation');
                    datatool_text = [datatool_text;chart_holder];
                end
                
                if ~isempty(datatool_inject)
                    % If there weren't any datacharts to begin with, append the DataCharts closer
                    datatool_text{end+1} = '</DataCharts>';
                end
                
                modified_text = [modified_text(1:datatool_inject0);...
                                datatool_text;...
                                modified_text(datatool_inject1:end)];
            end
            
            %%%Inject Stimulus Code
            try obj.stimulus_objects(1).name;            
                numStims = size(obj.stimulus_objects,1);
            catch
                numStims = 0;
            end 
            
            % Disable any existing stimuli
            extStims = find(contains(modified_text,'<Stimulus>'));
            temp_extStimHolder = {};
            if sum(contains(modified_text,'MotorPosition'))>0
                motorStimEnd = zeros(sum(contains(modified_text,'MotorPosition')),1);
            else
                motorStimEnd = [];
            end
            for ii = 1:length(extStims)
                stimEnd = find(contains(modified_text(extStims(ii):end),'</Stimulus>'),1,'first')+extStims(ii)-1;
                if ~isempty(find(contains(modified_text(extStims(ii):stimEnd),'MotorPosition')))
                    enInd = find(contains(modified_text(extStims(ii):end),'<Enabled>'),1,'first')+extStims(ii)-1;
                    modified_text{enInd} = '<Enabled>False</Enabled>';
                    motorStimEnd(ii) = stimEnd;
                    temp_extStimHolder = [temp_extStimHolder;modified_text(extStims(ii):stimEnd)];
                end
            end
            
            if ~isempty(extStims)
                modified_text = [modified_text(1:find(contains(modified_text,'<ExternalStimuli>')));...
                                temp_extStimHolder;...
                                modified_text(find(contains(modified_text,'</ExternalStimuli>')):end)];
            end
            
            if numStims > 0
%               stimuli_inject = find(contains(modified_text,'</ExternalStimuli>'));  
                if length(extStims) >= 1
                    if ~isempty(motorStimEnd)
                        stimuli_inject_0 = max(find(contains(modified_text,'</Stimulus>'),1,'last'));
                    else
                        stimuli_inject_0 = find(contains(modified_text,'<ExternalStimuli>')); 
                    end
                    stimuli_inject_1 = find(contains(modified_text,'</ExternalStimuli>')); 
                else
                    stimuli_inject_0 = find(contains(modified_text,'<ExternalStimuli/>'));
                    modified_text{stimuli_inject_0} = '<ExternalStimuli>';
                    stimuli_inject_1 = stimuli_inject_0+1;
                    modified_text = [modified_text(1:stimuli_inject_0);...
                                    '</ExternalStimuli>';...
                                    modified_text(stimuli_inject_1:end)];
                end
                stim_text = {};
                for i=1:numStims
                    stim_holder = CanvasText(obj).build_stimulus(obj.stimulus_objects(i),'simulation');
                    stim_text = [stim_text;stim_holder];
                end
%                 modified_text = [modified_text(1:stimuli_inject-1);...
%                                 stim_text;...
%                                 modified_text(stimuli_inject:end)];
                modified_text = [modified_text(1:stimuli_inject_0);...
                                stim_text;...
                                modified_text(stimuli_inject_1:end)];
            end
            
            sim_file_revised = strcat(sim_file(1:end-5),'_fake.asim');
            %sim_file_revised = [pwd,'\Animatlab\SynergyWalking\muscleStim.asim'];
            [~,simName,simExt] = fileparts(sim_file_revised);
            %disp(['Simulation file ',simName,simExt,' has been updated.'])
            fileID = fopen(sim_file_revised,'w');
            fprintf(fileID,'%s\n',modified_text{:});
            fclose(fileID);
        end
        %% update_link_cdata: Update the cdata numbers in the link objects
        function update_link_cdata(obj)
            no_links = 0;
            try no_links = isempty(fields(obj.link_objects));
            catch
                no_links  = 1;
            end
            if size(obj.neuron_objects,1) > 1 && ~no_links
                neurons = {obj.neuron_objects(:,1).name};
                for ii = 1:size(obj.link_objects,1)
                    link = obj.link_objects(ii);
                    origin_num = find(ismember(neurons, link.origin_ID(1:end-3)))-1;
                    destination_num = find(ismember(neurons, link.destination_ID(1:end-3)))-1;
                    obj.link_objects(ii).origin_cdata_num = origin_num;
                    obj.link_objects(ii).destination_cdata_num = destination_num;
                end
            end
        end
        %% ConstrainedPosition
        function pos = ConstrainedPosition(obj,type,pos,bounds)
            if strcmp(type,'n') %target
                sz = CanvasConstants.NEURON_size;
            elseif strcmp(type,'selectionbox')
                sz = 0;
            else %obstacle
                sz = CanvasConstants.OBSTACLE_SIZE;
            end
            
            bottomLeftLimits = [1 1] + sz/2;
            topRightLimits = bounds - [1 1] - sz/2;
            
            blCond = pos < bottomLeftLimits;
            pos(blCond) = bottomLeftLimits(blCond);
            
            urCond = pos > topRightLimits;
            pos(urCond) = topRightLimits(urCond);
            pos = floor(pos);
        end
        %% updateStimModel
        function updateStimModel(obj,neur_index)
                neuron = obj.neuron_objects(neur_index);
                dt = obj.dtsim;
                simendtime = obj.tmax;
                R = 20;
                Cm = .1;
                
                Upost = zeros(1,simendtime/dt+1);
                Upre = {};
                
                if ~isempty(neuron.inlinks)
                    Upre = cell(length(neuron.inlinks),3);
                    for i = 1:length(neuron.inlinks)
                        oNeur = obj.neuron_objects(contains({obj.neuron_objects.ID},neuron.inlinks(i).origin_ID));
                        oLink = neuron.inlinks(i);
                        oSyn = obj.synapse_types(contains({obj.synapse_types.name},oLink.synaptictype));
                        delE = oSyn.delE;
                        if ~isempty(oSyn.k) && ~isempty(oSyn.c)
                            keyboard
                            fprintf('Error: synapse has both a k and c value. Should have one or the other.')
                        else
                            if isempty(oSyn.k)
                                csyn = oSyn.c;
                                gMax = (csyn*R-R)/(delE-csyn*R);
                            else
                                ksyn = oSyn.k;
                                gMax = (ksyn*R)/(delE-ksyn*R);
                            end
                        end
                            Upre{i,1} = (dt/Cm)*(gMax/R);
                            Upre{i,2} = delE;
                        if ~isempty(oNeur.totmem)
                            Upre{i,3} = oNeur.totmem;
                        else
                            Upre{i,3} = zeros(1,simendtime/dt+1);
                        end
                    end
                end
                
                if ~isempty(neuron.stimulus)
                    %Hardcoded, to be changed
                    stim = neuron.stimulus;
                    dt = obj.dtsim;
                    simendtime = obj.tmax;

                    %Gather stimulus info
                    amp = stim.amplitude;
                    t = 0:dt:simendtime;

                    %Make new stimulus waveform
                    Iapp = zeros(size(t));
                    ind = @(time) time/dt+1;
                    start_ind = ind(stim.starttime);
                    end_ind = ind(stim.endtime);
                    Iapp(start_ind:end_ind) = amp;

                    %Store the stimulus waveform in the model
                    neuron.stimulus.waveform = Iapp;
                else
                    Iapp = zeros(1,simendtime/dt+1);
                end

                for i = 2:length(Upost)
                    inlink_adder = 0;
                    for j = 1:size(Upre,1)
                        inlink_adder = inlink_adder + Upre{j,1}*Upre{j,3}(i-1)*(Upre{j,2}-Upost(i-1));
                    end
                    stimterm(i-1) = (dt/Cm)*Iapp(i);
                    uposterm(i-1) = (1-dt/Cm)*Upost(i-1);
                    Upost(i) = (1-dt/Cm)*Upost(i-1)+(dt/Cm)*Iapp(i)+inlink_adder;
                end
                
                neuron.totmem = Upost;
                obj.neuron_objects(neur_index) = neuron;
        end
    end
end