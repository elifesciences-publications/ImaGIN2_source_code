function D = ImaGIN_spm_eeg_rdata_spike2_mono(S)
% Converts EEG data from CED Spike2 to SPM format (.mat/.dat)

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
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David


%% ===== INPUTS =====
% Prepare SPM window
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Finter)

% Parse inputs
try
    Fdata = S.Fdata;
catch
    Fdata = spm_select(1, '\.smr$', 'Select Spike2 file');
end
try
    FileOut = S.FileOut;
catch
    FileOut = spm_input('Name of converted file', '+1', 's');
end
try
    S.channel;
catch
    S.channel = str2num(spm_input('Index of channels to read', '+1', 's'));
end
try
    S.coarse;
catch
    S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
end
% Wait icon
spm('Pointer', 'Watch'); 
drawnow;


%% ===== READ FILE =====
% Read file header
fid = fopen(Fdata, 'r');
ChanList = SONChanList(fid);
fclose(fid);
% By default: read all channels
if isempty(S.channel)
    S.channel = [ChanList.number];
end
% Read data
[Data, iGoodChannels, Evt] = ImaGIN_read_smr(Fdata, S.channel, S.coarse);



%% ===== SAVE SPM FILE =====
% Separate path and file name
[fPath, fBase] = fileparts(FileOut);
% D structure
D.type      = 'continuous';
D.path      = fPath;
D.fname     = [fBase '.mat'];
D.Nchannels = size(Data.data,1);
D.Nsamples  = length(Data.time);
D.Fsample   = 1 ./ (Data.time(2) - Data.time(1));
D.time      = Data.time;

% Make channel names unique
chNames = Data.name;
for i = 1:length(chNames)
    if (nnz(strcmpi(chNames{i}, Data.name)) > 1)
        chNames{i} = sprintf('%s%02d', chNames{i}, i);
    end
end
% Channel names
for i = 1:length(chNames)
    % Channels
    D.channels(i).label = chNames{i};
    % Sensors
    D.sensors.eeg.label{i}       = D.channels(i).label;
    D.sensors.eeg.chantype{i}    = 'eeg';
    D.sensors.eeg.chanunit{i}    = 'm';
    D.sensors.eeg.chanpos(i,1:3) = [NaN, NaN, NaN];
    D.sensors.eeg.elecpos(i,1:3) = [NaN, NaN, NaN];
end
D.sensors.eeg.unit = 'm';


% Trials
D.trials.label  = 'Undefined';
D.trials.events = [];
D.trials.onset  = 1 ./ D.Fsample;

% Physically initialise .dat file
D.data = file_array([FileOut '.dat'], [D.Nchannels, D.Nsamples], 'float32-le');
D.data(:,:) = full(Data.data);

% Convert to SPM MEEG object
Dmeeg = meeg(D);
% Add events
if ~isempty(Evt)
    Dmeeg = events(Dmeeg, 1, Evt);
end
% Get full file
Dfile = fullfile(Dmeeg);

% Convert back to struct to be able to add the 'spike' field
D = struct(Dmeeg);
% Spikes
if isfield(Data,'spike')
    D.spike = Data.spike;
end
% Save file
save(Dfile, 'D');

% Stop processing
spm('Pointer','Arrow');


