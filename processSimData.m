function simStruct = processSimData(sim_path)
    % Convert a simulation file into a struct containing all data
    
    % Input: sim_path (char or string): absolute path for a simulation (.asim) file
    % Output: simStruct: struct containing all output data references in simulation file
    if nargin ~= 1
        error('Enter the path to an .asim file.')
    end

    if isstring(sim_path)
        sim_path = char(sim_path);
        if ~contains(sim_path,'.asim')
            error('Provided file path is not an .asim file.')
        end
    end
   
    if ~isfile(sim_path)
        error('Sim file doesn''t exist.')
    end
    
    % Pre-process .asim file to ensure that mesh filepaths and mass values are correct
%     meshMatch(sim_path);
%     massCheck(sim_path);
    
    % Run the simulation file
    nativeSimPath = 'C:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe';
    VSSimPath = 'C:\AnimatLabSDK\AnimatLabPublicSource\bin\AnimatSimulator.exe';
    whichSimulator2Use = [isfile(nativeSimPath) isfile(VSSimPath)];
    if whichSimulator2Use(2) == 1
        simulatorPath = VSSimPath;
    else
        if whichSimulator2Use(1) == 1
            simulatorPath = nativeSimPath;
        else
            error('No simulator detected. Is Animatlab installed on this system?')
        end
    end

    % Check is jsystem is on the Matlab file path. If not, add it.
    pathCell = regexp(path, pathsep, 'split');
    jSystemOnPath = any(strcmpi([pwd,'\Plugins\jsystem\src'], pathCell));
    
    if ~jSystemOnPath && isempty(ls('jsystem.m'))
        addpath([pwd,'\Plugins\jsystem\src'])
    end
    
    executable = [string(simulatorPath),string(sim_path)];
    [res,out,err] = jsystem(executable);
    if ~isempty(err)
        error(err)
    end
    
%     executable = ['"',simulatorPath,'" "',sim_path,'"'];
%     [status, message] = system(executable);
%     if status
%         error(message)
%         return
%     end

    % Process output simulation rsults, as stored in .txt files related to the .aform datatools
    sim_text = importdata(sim_path);
    temp_filenames = sim_text(contains(sim_text,'<OutputFilename>'));
    filenames = extractBetween(temp_filenames,'<OutputFilename>','</OutputFilename>');
    simStruct = struct();
    counter = 1;

    for ii = 1:size(filenames,1)
        fpath = [char(fileparts(sim_path)),'\',filenames{ii}];
        if isfile(fpath)
            ds = importdata(fpath);
            if ~isempty(ds)
                simStruct(counter).name = filenames{ii}(1:end-4);
                if isfield(ds,'colheaders')
                    simStruct(counter).time = ds.data(:,2);
                    simStruct(counter).data = ds.data(:,3:end);
                    simStruct(counter).data_cols = ds.colheaders(1,3:end);
                else
                    simStruct(counter).time = ds.textdata(2:end,2);
                    simStruct(counter).data = ds.data;
                end
                counter = counter+1;
            end
        else
            warning('%s does not exist.\n',filenames{ii})
        end
    end
end