file_dir = fileparts(mfilename('fullpath'));
proj_file = [pwd,'\Animatlab\SynergyWalking\SynergyWalking20200109.aproj'];
meshMatch(proj_file);
revised_file = strcat(proj_file(1:end-6),'.aproj');
[projDir,projName,ext] = fileparts(revised_file);
disp(['Starting to build Animatlab project ',projName,ext])
delete([projDir,'\Trace*'])

original_text = importdata(proj_file);
muscleIDs = find(contains(original_text,'<PartType>AnimatGUI.DataObjects.Physical.Bodies.LinearHillMuscle</PartType>'))-2;
muscle_out = scrape_project_for_musc_info(original_text);
numMuscles = length(muscleIDs);

nsys = CanvasModel;
neurpos = [];

nsys.addDatatool('MuscleLength');

for ii =1:38
    %scrape for muscle ii
    %add muscle to nsys
    nsys.addMuscle([muscle_out{ii,2},'-neural'],muscle_out{ii,3},[ii+10 50])
    %add DTaxes for that muscle
    nsys.addDTaxes(nsys.datatool_objects(1),nsys.muscle_objects(ii),'MuscleLength')
end

nsys.create_animatlab_project(proj_file);
disp(['Animatlab project file ',projName,ext,' created.'])

function equation = generate_synergy_eq(bigH,simTime)
    dt = .00054;
%     time = (99*dt:dt:10.01-10*dt)';
    time = (0:dt:simTime)';
    
    bigH = interpolate_for_time(time,bigH);
    
    % Create a sum of sines equation for the joint angle waveforms
    fitresult = sumsinesFit(time, bigH,8);
    % Coeffs are the a, b, and c values in the equation a*sin(b*t+c)
    coeffs = [coeffnames(fitresult),num2cell(coeffvalues(fitresult)')];
    % Equations are in the format necessary for integration into Animatlab's .asim filetype
    equation = sum_of_sines_maker(coeffs,1);
end

function waveformsBig = interpolate_for_time(time,waveforms)
    % Interpolate the undersampled input to match the required time vector
    waveforms = waveforms';
    m = length(time);
    n = length(waveforms);
    if m ~= n
        waveformsBig = interp1(1:n,waveforms,linspace(1,n,m));
    end

    avgblocks = floor(.01*length(waveformsBig));
    coeffblocks = ones(1,avgblocks)/avgblocks;
    for i=1:2
        waveformsBig = filtfilt(coeffblocks,1,waveformsBig);
    end
    waveformsBig = waveformsBig';
end

function muscle_out = scrape_project_for_musc_info(input_text)
    % Scrape text for muscle information
    % Input: input_text: cell array of text from an Animatlab project
    % Output: muscle_out: cell array containing muscle indices, muscle names, and muscle ID's
    muscleInds = find(contains(input_text,'<PartType>AnimatGUI.DataObjects.Physical.Bodies.LinearHillMuscle</PartType>'))-2;
    if isempty(muscleInds)
        error('No muscles in input text')
    end
    muscle_out = cell(length(muscleInds),3);
    for ii = 1:length(muscleInds)
        muscle_out{ii,1} = muscleInds(ii);
        muscle_out{ii,2} = lower(input_text{muscleInds(ii)-1}(7:end-7));
        muscle_out{ii,3} = input_text{muscleInds(ii)}(5:end-5);
    end
end