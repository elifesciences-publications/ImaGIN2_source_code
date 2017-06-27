function tutorial_epimap()
% TUTORIAL_EPIMAP Script corresponding to the ImaGIN/epileptogenicity tutorial.
%
% DESCRIPTION:
%     Corresponding online tutorial: [URL]
%
%     This example script processes only one subject, but illustrates the 
%     structure corresponding to a study with multiple subjects.

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2017 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David,  2010-2017
%          Francois Tadel, 2017


%% ===== DATA DEFINITION =====
% Define tutorial folder
% Root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
Root = 'C:\Work\Dev\Inserm';

% Patient #1: Input files and processing options
I = 1;
% Name of the patient, which must corresponds to subfolder in PatientsDir
Patient{I}.Name = 'tutorial_epimap';
% Patient anatomy: pre- and post-implantation MRI scans
Patient{I}.MRI.pre  = fullfile(Root, Patient{I}.Name, 'anat', 'MRI', '3DT1pre_deface.nii');
Patient{I}.MRI.post = fullfile(Root, Patient{I}.Name, 'anat', 'MRI', '3DT1post_deface.nii');
Patient{I}.MRI.ref  = 'pre';
Patient{I}.MRI.out  = fullfile(Root, Patient{I}.Name, 'anat', 'MRI');
% Name of the files with the seizure recordings (original Micromed have a .TRC extension)
Patient{I}.File{1}     = 'SZ1';   % Short seizure, no propagation   
Patient{I}.File{2}     = 'SZ2';   % Short seizure, propagation 
Patient{I}.File{3}     = 'SZ3';   % Long seizure, generalized 
% Seizure onset, from the beginning of the file (the events "Seizure" indicated in the .TRC files are not reliable)
Patient{I}.Onset = [120.800, ...         % File #1:  423ms after the Seizure marker
                    143.510, ...         % File #2: 2900ms after the Seizure marker
                    120.287];            % File #3:   60ms after the Seizure marker
                
% % Baseline segment (20s), with respect to the Onset marker: Defined while reviewing the recordings as a bipolar montage
% Patient{I}.Baseline = {[-48, -28], ...   % File #1: From beginning of recordings 72.8s-92.8s
%                        [-40, -20], ...   % File #2: From beginning of recordings 103.5s-123.5s
%                        [-75, -55]};      % File #3: From beginning of recordings 45.3s-65.3s
% Shorter baseline for faster execution: 5s
Patient{I}.Baseline = {[-48, -43], ...   % File #1: From beginning of recordings 72.8s-92.8s
                       [-40, -35], ...   % File #2: From beginning of recordings 103.5s-123.5s
                       [-75, -70]};      % File #3: From beginning of recordings 45.3s-65.3s
Patient{I}.BaselineFile = {[],[],[]};

% List of bad channels: Defined while reviewing the recordings as a bipolar montage
Patient{I}.BadChannel  = {[74 28], ...   % File #1: v'2v'1, f'2f'1
                          [74], ...      % File #2: v'2v'1
                          [74 54]};      % File #3: v'2v'1, o'2o'1
% Epileptogenicity options
Patient{I}.FreqBand     = [120 200];  % Defined by looking at the TF maps (using the same for the three seizures)
Patient{I}.TimeConstant = 3;          % Duration of the sliding window of interest: 3s
Patient{I}.Latency      = 0:2:20;     % For files #2 and #3, will compute delay maps with sliding windows of 3s between 0s and 20s post-seizure
Patient{I}.Prefix       = '';
ThDelay = 0.05;


%% ===== PREPARE ANATOMY =====
% Co-registration, segmentation and normalization of the pre-op and post-op MRI scans
ImaGIN_anat_spm(Patient);


%% ===== IMPORT SEEG =====
% Loop on subjects (only on this example)
for i0 = 1:length(Patient)
    % Delete previously created files
    fileMat = dir(fullfile(Root, Patient{i0}.Name, 'seeg', '*.mat'));
    fileDat = dir(fullfile(Root, Patient{i0}.Name, 'seeg', '*.dat'));
    fileTxt = dir(fullfile(Root, Patient{i0}.Name, 'seeg', '*.txt'));
    if ~isempty(fileMat) || ~isempty(fileDat) || ~isempty(fileTxt)
        AllFiles = cellfun(@(c)fullfile(Root, Patient{i0}.Name, 'seeg', c), {fileMat.name, fileDat.name, fileTxt.name}, 'UniformOutput', 0);
        delete(AllFiles{:});
    end
    
    % Loop on seizure datasets
    for i1 = 1:length(Patient{i0}.File)
        % Convert Micromed .TRC to SPM .mat/.dat
        clear S
        S.dataset = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.TRC']);
        S.FileOut = fullfile(Root, Patient{i0}.Name, 'seeg', Patient{i0}.File{i1});
        D = ImaGIN_spm_eeg_converteeg2mat(S);

        % Add electrodes positions
        clear S
        S.Fname        = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        S.filenamePos  = fullfile(Root, Patient{i0}.Name, 'anat', 'implantation', 'Electrodes_Pos_MNI.txt');
        S.filenameName = fullfile(Root, Patient{i0}.Name, 'anat', 'implantation', 'Electrodes_Name.txt');
        D = ImaGIN_Electrode(S);
        
        % Longitudinal bipolar montage
        clear S
        S.Fname    = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        S.FileOut  = fullfile(Root, Patient{i0}.Name, 'seeg', ['b' Patient{i0}.File{i1} '.mat']);
        D = ImaGIN_BipolarMontage(S);
        [tmp, Patient{i0}.FileBipolar{i1}] = fileparts(S.FileOut);
    end
end


%% ===== REVIEW =====
% Review the recordings for each subject and each seizure:
%   - Mark the onset of the seizure (in these files, the start of the seizures is already marked in the original TRC files)
%   - Identify bad channels.

% % Example file and channels to display
% FilaName = fullfile(Root, Patient{1}.Name, 'seeg', [Patient{1}.FileBipolar{1} '.mat']);

% % ImaGIN viewer
% SelChan = 1:5;
% ImaGIN_DispData(FilaName, SelChan);

% % SPM viewer
% D = spm_eeg_load(FilaName);
% spm_eeg_review(D);

% % Display contact positions
% S.Fname = fullfile(Root, Patient{1}.Name, 'seeg', [Patient{1}.FileBipolar{1} '.mat']);
% S.P     = fullfile(Root, Patient{1}.Name, 'anat', 'MRI', 'wBrainPre.nii');
% ImaGIN_DispElectrodes(S);




%% ===== AFTER REVIEW =====
for i0 = 1:length(Patient)
    for i1 = 1:length(Patient{i0}.FileBipolar)
        FileName = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{i1} '.mat']);

        % ===== SET BAD CHANNELS =====
        if ~isempty(Patient{i0}.BadChannel)
            % % SPM version
            % D = spm_eeg_load(FileName);
            % D = badchannels(D, Patient{i0}.BadChannel{i1}, 1);
            % save(D);
            % ImaGIN version
            clear S
            S.Fname       = FileName;
            S.BadChannels = Patient{i0}.BadChannel{i1};
            ImaGIN_BadChannelSet(S);
        end
        
        % ===== SET ONSET EVENTS =====
        clear S
        S.Fname        = FileName;
        S.Action       = 'Add';
        S.Nevent       = 1;
        S.EventName{1} = 'Onset';
        S.Timing{1}    = Patient{I}.Onset(i1);
        ImaGIN_Events(S);

        % ===== SET TIME ORIGIN =====
        clear S
        S.Fname    = FileName;
        S.EventRef = 'Onset';
        S.Offset   = 0;
        ImaGIN_TimeZero(S);
        
        % ===== IMPORT BASELINE ======
        if ~isfield(Patient{i0}, 'BaselineFile') || (length(Patient{i0}.BaselineFile) < i1) || isempty(Patient{i0}.BaselineFile{i1})
            clear S
            S.Fname       = FileName;
            S.Job         = 'Manual';
            S.EventStart  = 'Onset';
            S.EventEnd    = 'Onset';
            S.OffsetStart = -Patient{i0}.Baseline{i1}(1);
            S.OffsetEnd   = Patient{i0}.Baseline{i1}(2);
            S.NewFile     = 1;
            S.Prefix      = 'Baseline_';
            ImaGIN_Crop(S);
            % Save name of the output file
            Patient{i0}.BaselineFile{i1} = [S.Prefix, Patient{i0}.FileBipolar{i1}];
        end
    end
end


%% ===== TIME-FREQUENCY: MULTITAPER =====
% Averaging the three seizures together
for i0 = 1:length(Patient)
    % Get time window to process (at most 10s before and after t=0)
    minTime = -10;
    maxTime = 10;
    for i1 = 1:length(Patient{i0}.FileBipolar)
        clear SS
        SS.D = fullfile(Root,Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{i1} '.mat']);
        D = spm_eeg_load(SS.D);
        minTime = max([minTime, min(time(D)) + 0.5]);
        maxTime = min([maxTime, max(time(D)) - 0.5]);
    end
    
    % Compute the TF decomposition with multi-taper
    for i1 = 1:length(Patient{i0}.FileBipolar)
        % Multi-taper decomposition
        clear SS
        SS.D = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{i1} '.mat']);
        D = spm_eeg_load(SS.D);
        SS.Pre             = '';
        SS.Method          = 'Multitaper';
        SS.frequencies     = 10:3:230;
        SS.FactMod         = 10;
        SS.TimeWindowWidth = 1;
        SS.TimeWindow      = [minTime maxTime];
        SS.TimeResolution  = 0.1;
        SS.NSegments       = 1;
        SS.Taper           = 'Hanning';
        SS.channels        = 1:D.nchannels;
        ImaGIN_spm_eeg_tf(SS);
        % Baseline normalization of the TF maps
        clear SS2
        SS2.D = fullfile(D.path,['m1_' SS.Pre '_' D.fname]);
        SS2.B = [-10 -1];
        ImaGIN_NormaliseTF(SS2);
    end
    % Average the TF maps
    [files,dirs] = spm_select('List', fullfile(Root, Patient{i0}.Name, 'seeg'), '^nm1.*\.mat$');
    if (size(files,1) > 1)
        clear S;
        S.D       = [repmat([fullfile(Root, Patient{i0}.Name, 'seeg'), filesep], size(files,1), 1), files];
        S.Method  = 'Mean';
        S.NewName = 'mean_tf_multitaper';
        D = ImaGIN_AverageTF(S);
        % Display results
        % ImaGIN_DispTF(D);
        % ImaGIN_DispTF(fullfile(Root, Patient{i0}.Name, 'seeg', [S.NewName, '.mat']));
    end
end


%% ===== EPILEPTOGENICITY: SEIZURE #1 =====
% Process separately seizure #1 (no propagation) and seizures #2 and #3 (generalized)
i0 = 1;
clear S;
% List of input files
S.D = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{1} '.mat']);            % Seizure data
S.B = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.BaselineFile{1} '.mat']);    % Baseline data
% Process options
S.TimeWindow     = (0 : 0.01 : Patient{i0}.TimeConstant+1+max(Patient{i0}.Latency));
S.FreqBand       = Patient{i0}.FreqBand;
S.HorizonT       = Patient{i0}.TimeConstant;
S.BadChannel     = Patient{i0}.BadChannel{1};
S.Latency        = 0;        % No propagation: Study only the first 3s (TimeConstant) after the Onset marker
S.TimeResolution = 0.2;
S.ThDelay        = ThDelay;
S.Atlas          = 'Human';
S.AR             = 0;
S.sMRI           = fullfile(Patient{I}.MRI.out, 'wBrainPre.nii');
S.CorticalMesh   = 1;
S.FileName       = Patient{i0}.Prefix;
ImaGIN_Epileptogenicity(S);

% Display of epileptogenicity map


    
%% ===== EPILEPTOGENICITY: SEIZURE #2-3 =====
clear S;
% List of input files
S.D = char(fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{2} '.mat']), ...          % Seizure data
           fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.FileBipolar{3} '.mat'])); 
S.B = char(fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.BaselineFile{2} '.mat']), ...  % Baseline data
           fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.BaselineFile{3} '.mat']));
% Process options
S.TimeWindow     = (0 : 0.01 : Patient{i0}.TimeConstant+1+max(Patient{i0}.Latency));
S.FreqBand       = Patient{i0}.FreqBand;
S.HorizonT       = Patient{i0}.TimeConstant;
S.BadChannel     = [Patient{i0}.BadChannel{2:3}];
S.Latency        = Patient{i0}.Latency;        % Use the sliding windows defined for this subject
S.TimeResolution = 0.2;
S.ThDelay        = ThDelay;
S.Atlas          = 'Human';
S.AR             = 0;
S.sMRI           = fullfile(Patient{I}.MRI.out, 'wBrainPre.nii');
S.CorticalMesh   = 1;
S.FileName       = Patient{i0}.Prefix;
ImaGIN_Epileptogenicity(S);



