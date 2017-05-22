function tutorial_epimap()
% TUTORIAL_EPIMAP Script correponding to the ImaGIN/epileptogenicity tutorial.
%
% DESCRIPTION:
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
% Get tutorial folder
Root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
% Patient #1: Input files and processing options
I = 1;
Patient{I}.Name        = 'tutorial_epimap';  % Name of the subject, which must corresponds to subfolder in PatientsDir
Patient{I}.File{1}     = 'SZ1';              % Name of the files with the seizure recordings (original Micromed have a .TRC extension)
Patient{I}.File{2}     = 'SZ2';
Patient{I}.File{3}     = 'SZ3';
Patient{I}.Baseline    = {};
% Patient{I}.Baseline{1} = 'Baseline_SZ1';   % Name of the files for the corresponding baseline recordings
% Patient{I}.Baseline{2} = 'Baseline_SZ2';
% Patient{I}.Baseline{3} = 'Baseline_SZ3';
Patient{I}.BadChannel  = {[74], ...          % List of bad channels  - File #1 => v'1
                          [74 28], ...         %  - File #2  => v'1, f'1
                          [74 54]};            %  - File #3  => o'1
% ====== TODO: weird to define bad channels after bipolar montage no?? =====
Patient{I}.FreqBand    = [210 230];          % TODO: Or [100 180] ???
Patient{I}.Latency     = 0;
Patient{I}.TimeConstant= 2;
Patient{I}.sMRI        = fullfile(Root, Patient{I}.Name, 'anat', 'processed', 'wBrainPre.nii');
Patient{I}.Pre         = '';
% Epileptogenicity options
ThDelay = 0.05;


%% ===== IMPORT =====
% Loop on subjects (only on this example)
for i0 = 1:length(Patient)
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
        S.SaveFile = Patient{i0}.File{i1};
        S.FileOut  = S.Fname;
        D = ImaGIN_BipolarMontage(S);
    end    
end


%% ===== REVIEW =====
% Review the recordings for each subject and each seizure:
%   - Mark the onset of the seizure (in these files, the start of the seizures is already marked in the original TRC files)
%   - Identify bad channels.


%% ===== AFTER REVIEW =====
for i0 = 1:length(Patient)
    for i1 = 1:length(Patient{i0}.File)
        % ===== SET TIME ORIGIN =====
        clear S
        S.Fname    = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        S.EventRef = 'Seizure';
        S.Offset   = 0;
        ImaGIN_TimeZero(S);
        
        % ===== SET BAD CHANNELS =====
        if ~isempty(Patient{i0}.BadChannel)
            D = spm_eeg_load(fullfile(Root,Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']));
            D = badchannels(D, Patient{i0}.BadChannel{i1}, 1);
            save(D);
        end
        
        % ===== IMPORT BASELINE ======
        clear S
        S.Fname       = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        S.Job         = 'Manual';
        S.EventStart  = [];
        S.EventEnd    = [];
        S.OffsetStart = 0;
        S.OffsetEnd   = 20;
        S.NewFile     = 1;
        S.Prefix      = 'Baseline_';
        ImaGIN_Crop(S);
        % Save name of the output file
        Patient{i0}.Baseline{i1} = [S.Prefix, Patient{i0}.File{i1}];
    end
end


% %% ===== TIME-FREQUENCY: WAVELET =====
% for i0 = 1:length(Patient)
%     for i1=1:length(Patient{i0}.File)
%         % Wavelet decomposition
%         clear SS
%         SS.D = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
%         D = spm_eeg_load(SS.D);
%         SS.Synchro         = 'No';
%         SS.Pre             = '';
%         SS.Method          = 'Morlet wavelet';
%         SS.frequencies     = 10:3:230;
%         SS.FactMod         = 0;
%         SS.Mfactor         = 20;
%         SS.Width           = 0;
%         SS.TimeWindow      = -10:.2:10;
%         SS.TimeWindowWidth = SS.TimeWindow(2) - SS.TimeWindow(1);
%         SS.Coarse          = 0;
%         SS.channels        = 1:D.nchannels;
%         SS.TimeResolution  = 0.05;
%         ImaGIN_spm_eeg_tf(SS);
%         % Baseline normalization of the TF maps
%         clear SS2
%         SS2.D=fullfile(D.path,['w1_' SS.Pre '_' D.fname]);
%         SS2.B=[-10 -1];
%         ImaGIN_NormaliseTF(SS2);
%     end
%     % Average the TF maps
%     [files,dirs] = spm_select('List', fullfile(Root, Patient{i0}.Name, 'seeg'), '^nw1.*\.mat$');
%     if (size(files,1) > 1)
%         clear S;
%         S.D       = [repmat([fullfile(Root, Patient{i0}.Name, 'seeg'), filesep], size(files,1), 1), files];
%         S.Method  = 'Mean';
%         S.NewName = 'mean_tf_wavelet';
%         D = ImaGIN_AverageTF(S);
%     end
% end


%% ===== TIME-FREQUENCY: MULTITAPER =====
for i0 = 1:length(Patient)
    % Get time window to process (at most 10s before and after t=0)
    minTime = -10;
    maxTime = 10;
    for i1 = 1:length(Patient{i0}.File)
        clear SS
        SS.D = fullfile(Root,Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        D = spm_eeg_load(SS.D);
        minTime = max([minTime, min(time(D)) + 0.5]);
        maxTime = min([maxTime, max(time(D)) - 0.5]);
    end
    
    % Compute the TF decomposition with multi-taper
    for i1 = 1:length(Patient{i0}.File)
        % Multi-taper decomposition
        clear SS
        SS.D = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
        D = spm_eeg_load(SS.D);
        SS.Pre             = '';
        SS.Method          = 'Multitaper';
        SS.Taper           = 'Hanning';
        SS.TimeResolution  = 0.1;
        SS.frequencies     = 10:3:230;
        SS.FactMod         = 10;
        SS.TimeWindow      = [minTime maxTime];
        SS.TimeWindowWidth = 1;
        SS.channels        = 1:D.nchannels;
        SS.NSegments       = 1;
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
    end
end


%% ===== EPILEPTOGENICITY =====
for i0 = 1:length(Patient)
    % Prepare list of input files
    clear S
    for i1 = 1:length(Patient{i0}.File)
        if i1 ==1
            S.D = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']);
            S.B = fullfile(Root, Patient{i0}.Name, 'seeg', [Patient{i0}.Baseline{i1} '.mat']);
        else
            S.D = char(S.D, fullfile(Root,Patient{i0}.Name, 'seeg', [Patient{i0}.File{i1} '.mat']));
            S.B = char(S.B, fullfile(Root,Patient{i0}.Name, 'seeg', [Patient{i0}.Baseline{i1} '.mat']));
        end
    end
    S.TimeWindow = (0 : 0.01 : Patient{i0}.TimeConstant+1+max(Patient{i0}.Latency));
    S.FreqBand   = Patient{i0}.FreqBand;
    S.HorizonT   = Patient{i0}.TimeConstant;
    S.BadChannel = Patient{i0}.BadChannel{i1};
    S.Latency        = Patient{i0}.Latency;
    S.TimeResolution = 0.1;
    S.ThDelay        = ThDelay;
    S.Atlas          = 'Human';
    S.AR             = 0;
    S.Latency        = 0;
    S.sMRI           = Patient{i0}.sMRI;
    S.CorticalMesh   = 1;
    try
        S.FileName = Patient{i0}.Pre;
    catch
        S.FileName = '';
    end
    ImaGIN_Epileptogenicity(S);
end


