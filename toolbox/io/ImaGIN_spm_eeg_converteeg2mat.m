function D = ImaGIN_spm_eeg_converteeg2mat(S)
% User interface for conversion of EEG-files to ImaGIN/SPM's data structure
%
% USAGE:  D = spm_eeg_converteeg2mat(S)
%         D = spm_eeg_converteeg2mat()
%
% INPUTS: 
%    - S: Optional struct with the followint (optional) fields:
%         |- dataset        : Full path of the file to convert
%         |- SelectChannels : Channels to keep in the conversion:
%         |                   Array of indices ([1,2,...]), Array of cells ({'A1','A2',...}) or String ('A1 A2 ...')
%         |                   If empty: the list of selected channels is auto-detected based on standard setups
%         |- FileOut        : Full path to the output .mat/.dat file (without file extension)
%         |- isSEEG         : 0 if importing surface EEG, 1 if importing SEEG/ECOG 
%         |- Fdata
%         |- Fchannels - String containing name of channel template file
%         | ... other format-specific fields

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
% Authors: Stefan Kiebel,  2005
%          Olivier David,  2010-2017
%          Francois Tadel, 2017

% Parse inputs
S2 = struct();
if (nargin < 1) || isempty(S)
    S = struct();
else
    S2 = ImaGIN_copy_fields(S2, S, {'dataset', 'Atlas', 'FileOut'});
end

% Ask dataset (if not defined in input)
if ~isfield(S2, 'dataset')
    [S.dataset, sts] = spm_select(1, '.*', 'Select M/EEG data file');
    S2.dataset = S.dataset;
    if ~sts
        return;
    end
end
% Selected channels: if not specified, auto-detect
if isfield(S, 'SelectChannels')
    SelectChannels = S.SelectChannels;
else
    SelectChannels = str2num(spm_input('Index of channels (empty=detect)', '+1', 's'));
end
% Ask if the input data is SEEG
if isfield(S, 'isSEEG') && ~isempty(S.isSEEG)
    isSEEG = S.isSEEG;
else
    isSEEG = isequal(spm_input('SEEG? ','+1','Yes|No'), 'Yes');
end

% Is output file defined
isOutputSet = isfield(S2, 'FileOut') && ~isempty(S2.FileOut);

% Processing starts
spm('Pointer','Watch'); drawnow;

% Loop on input files
Nfiles = size(S.dataset, 1);
D = cell(1,Nfiles);
for i1 = 1:Nfiles
    % Input file
    S2.Fdata = deblank(S.dataset(i1, :));
    % Get extension
    [fPath, fBase, fExt] = fileparts(S2.Fdata);
    % Default output filename
    if ~isOutputSet
        S2.FileOut = fullfile(fPath, fBase);
    end
    
    % Old versions of the readers
    % BrainAmp .eeg
        % S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'Bipolar', 'Montage', 'epochlength', 'coarse', 'channel', 'SaveFile'}); % MISSING FUNCTION EEG2MAT
        % D{i1} = ImaGIN_spm_eeg_rdata_elan(S2);
    % Micromed .trc
        % Old version: JPL & OD
        % S2 = ImaGIN_copy_fields(S2, S, {'pts', 'System', 'CreateTemplate', 'Fchannels', 'channel', 'coarse', 'Montage', 'MontageName', 'NeventType', 'event_file' , 'Bipolar', 'bipole', 'loadevents'});
        % D{i1} = ImaGIN_spm_eeg_rdata_micromed_mono(S2);
    % EDF
        % S2 = ImaGIN_copy_fields(S2, S, {'channel', 'Atlas', 'SEEG', 'Bipolar', 'coarse', 'SizeMax'});
        % D = ImaGIN_spm_eeg_rdata_edf(S2);
    % Nicolet .e
        % S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Montage', 'coarse', 'SEEG', 'filenamePos', 'filenameName', 'MontageName', 'SaveFile', 'channel'});
        % D = ImaGIN_spm_eeg_rdata_nicolet_mono(S2);
    % Deltamed .bin
        % S2 = ImaGIN_copy_fields(S2, S, {'Bipolar', 'Bipole', 'CreateTemplate', 'Montage', 'coarse', 'Nevent', 'SEEG', 'filenamePos', 'filenameName', 'MontageName', 'SaveFile'});
        % D{i1} = ImaGIN_spm_eeg_rdata_deltamedbin_mono(S2);
            
    % Switch between file formats
    BstFormat = [];
    switch lower(fExt)
        case '.smr'
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'Bipolar', 'Montage', 'epochlength', 'coarse', 'channel'});
            D{i1} = ImaGIN_spm_eeg_rdata_spike2_mono(S2);
            
        case {'.asc','.txt'}
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'coarse', 'Bipolar', 'channel', 'Radc', 'Nevent', 'filenamePos', 'filenameName', 'MontageName'});
            D{i1} = ImaGIN_spm_eeg_rdata_ascii(S2);

        case {'.msm', '.msr'}
            S2 = ImaGIN_copy_fields(S2, S, {'Atlas', 'SEEG', 'Bipolar', 'coarse', 'SaveFile'});
            D{i1} = ImaGIN_spm_eeg_rdata_msm_mono(S2);

        case '.eeg'  % BrainAmp or Nihon Kohden
            % BrainAmp: There is a header in the same folder (.vhdr or .ahdr)
            if exist(fullfile(fPath, [fBase, '.vhdr']), 'file') || exist(fullfile(fPath, [fBase, '.ahdr']), 'file')
                BstFormat = 'EEG-BRAINAMP';
            % Nihon Kohden
            else
                BstFormat = 'EEG-NK';
            end
            
        case '.bin'
            BstFormat = 'EEG-DELTAMED';
        case '.trc'
            BstFormat = 'EEG-MICROMED';
        case '.edf'
            BstFormat = 'EEG-EDF';
        case '.e'
            BstFormat = 'EEG-NICOLET';
        otherwise
            error('Unknown format');
    end
    
    % Read file with Brainstorm functions
    if ~isempty(BstFormat)
        D{i1} = ImaGIN_convert_brainstorm(S2.Fdata, BstFormat, [S2.FileOut, '.mat'], SelectChannels, isSEEG);
    end
end
% Processing stops
spm('Pointer','Arrow');
end



%% ===== CONVERT WITH BRAINSTORM I/O =====
% Converts EEG data from a native file to SPM-format, and saves a log of the channel selection
function D = ImaGIN_convert_brainstorm(InputFile, FileFormat, OutputFile, SelChannels, isSEEG) %#ok<STOUT>
    % Open file with Brainstorm fuctions
    switch (FileFormat)
        % MICROMED .TRC
        case 'EEG-MICROMED'
            [sFileIn, ChannelMat] = in_fopen_micromed(InputFile);
            % Fix the names of some channels
            % Code copied from original function: ImaGIN_spm_eeg_rdata_micromed_mono
            for i = 1:length(ChannelMat.Channel)
                chName = ChannelMat.Channel(i).Name;
                % Lyon
                chName = lower(chName(~isspace(chName)));
                % Salpetriere: Can have both "C 4" and "c4" in the file: If a channel is duplicated, add a "s" after (the scalp EEG is always at the end)
                if ismember(chName, {ChannelMat.Channel(1:i-1).Name})
                    chName = [chName, 's'];
                end
                % Rennes
                chName = chName(setdiff(1:length(chName), strfind(chName,'.')));
                % Remove zeros
                v_num = find((chName=='0') | (chName=='1') | (chName=='2') | (chName=='3') | (chName=='4') | (chName=='5') | (chName=='6') | (chName=='7') | (chName=='8') | (chName=='9'));
                if (length(v_num) > 1)
                    if (v_num(2) == v_num(1) + 1) && strcmp(chName(v_num(1)), '0')
                        chName = chName(setdiff(1:length(chName), v_num(1)));
                    end
                end
                ChannelMat.Channel(i).Name = chName;
            end
        % DELTAMED .bin
        case 'EEG-DELTAMED'
            [sFileIn, ChannelMat] = in_fopen_deltamed(InputFile);
        % NIHON KOHDEN .EEG
        case 'EEG-NK'
            [sFileIn, ChannelMat] = in_fopen_nk(InputFile);
        % EDF / XLTEK
        case 'EEG-EDF'
            [sFileIn, ChannelMat] = in_fopen_edf(InputFile);
            
            % For XLTEK files exported as EDF + events in .txt file: reading the events
            XltekEvtFile = [InputFile(1:end-3) 'txt'];
            if exist(XltekEvtFile, 'file')
                % Read events
                sFileIn.events = in_events_xltek(sFileIn, XltekEvtFile);
                % If there are more than 3 events: do some cleaning
                if (length(sFileIn.events) > 3)
                    % Remove the first 3 event types
                    sFileIn.events(1:3) = [];
                    % Remove some other event types based on their names
                    iDel = find(ismember({sFileIn.events.label}, {'XLEvent','XLSpike','Gain/Changement de filtre','e','h','s','Opened relay','Closed relay'}) | ...
                                ~cellfun(@(c)isempty(strfind(c, 'Closed relay')), {sFileIn.events.label}));
                    if ~isempty(iDel)
                        sFileIn.events(iDel) = [];
                    end
                end
            end
            
            % For recordings from YUQ: many channels are classified as "POL" => Changing the type manually to EEG
            iPOL = find(strcmpi({ChannelMat.Channel.Type}, 'POL'));
            if ~isempty(iPOL)
                [ChannelMat.Channel(iPOL).Type] = deal('EEG');
            end
            
        % BrainVision BrainAmp
        case 'EEG-BRAINAMP'
            [sFileIn, ChannelMat] = in_fopen_brainamp(InputFile);
        % Nicolet
        case 'EEG-NICOLET'
            [sFileIn, ChannelMat] = in_fopen_nicolet(InputFile);
        otherwise
            error(['Unsupported file format: ', FileFormat]);
    end
    
    % Get list of available channel names
    ChanLabelsIn = {ChannelMat.Channel.Name};
    % Auto-detect good SEEG channels
    if isempty(SelChannels)
        % Get channels classified as EEG
        iEEG = channel_find(ChannelMat.Channel, 'EEG,SEEG,ECOG,ECG,EKG');
        % If there are no channels classified at EEG, take all the channels
        if isempty(iEEG)
            disp('ImaGIN> Warning: Channel types are not defined: not selecting by types.');
            iEEG = 1:length(ChannelMat.Channel);
        end
        % Detect channels of interest
        [iSelEeg, iEcg] = ImaGIN_select_channels({ChannelMat.Channel(iEEG).Name}, isSEEG);
        % If no SEEG channels, read as EEG
        if isempty(iSelEeg) && (isSEEG == 1)
            disp('ImaGIN> No SEEG channels selected by name, reading as regular EEG instead.');
            iSelEeg = ImaGIN_select_channels({ChannelMat.Channel(iEEG).Name}, 0);
            isSEEG = 0;
        end
        % Convert indices back to the original list of channels
        if isempty(iSelEeg)
            disp('ImaGIN> No channels selected by name, reading as regular EEG instead of SEEG.');
            iSel = iEEG;
            isSEEG = 0;
        else
            iSel = iEEG(iSelEeg);
        end
    % Selected channels are passed in input
    else
        % Find channels: string, cell array of strings, or array of indices
        iSel = ImaGIN_find_channels(ChanLabelsIn, SelChannels);
        iEcg = [];
    end
    % If no channels are selected
    if isempty(iSel)
        error('No valid channel names were found.');
    end
    % Set channels to ECG
    if ~isempty(iEcg)
        [ChannelMat.Channel(iEcg).Type] = deal('ECG');
    end
    % Keep only these ones in the data
    ChannelMat.Channel = ChannelMat.Channel(iSel);
   
    % Get number of epochs
    if isempty(sFileIn.epochs)
        nEpochs = 1;
    else
        nEpochs = length(sFileIn.epochs);
    end
    % Loop on the epochs
    for iEpoch = 1:nEpochs
        % Read entire segment with Brainstorm functions
        F = in_fread(sFileIn, ChannelMat, iEpoch, []);
        F = F(iSel,:);

        % File pointer to the output SPM file: keep only on epoch and a subset of the channels
        sFileInSpm = sFileIn;
        sFileInSpm.channelflag  = sFileIn.channelflag(iSel);
        % Remove all the events that are not in this segment
        if (nEpochs > 1)
            sFileInSpm.prop.samples = sFileIn.epochs(iEpoch).samples;
            sFileInSpm.prop.times   = sFileIn.epochs(iEpoch).times;
            for iEvt = 1:length(sFileInSpm.events)
                iEvtEpoch = find(sFileInSpm.events(iEvt).epochs == iEpoch);
                sFileInSpm.events(iEvt).times   = sFileInSpm.events(iEvt).times(iEvtEpoch);
                sFileInSpm.events(iEvt).samples = sFileInSpm.events(iEvt).samples(iEvtEpoch);
                sFileInSpm.events(iEvt).epochs  = sFileInSpm.events(iEvt).epochs(iEvtEpoch);
            end
            sFileInSpm.epochs = [];
        end
        
        % Output name of the SPM .mat/.dat file
        if (nEpochs == 1)
            SpmFile = OutputFile;
        else
            SpmFile = strrep(OutputFile, '.mat', sprintf('-%d.mat', iEpoch));
        end
        % Export to SPM format
        sFileOut = out_fopen_spm(SpmFile, sFileInSpm, ChannelMat);
        out_fwrite_spm(sFileOut, [], [], F);
        % Load new file to return the D structure
        load(SpmFile);

        % Save the original list of channels in a log file
        ImaGIN_save_log(SpmFile, ['Convert: Available channels     (' InputFile ')'], ChanLabelsIn);
        % Save the list of channels in the output file
        ChanLabelsOut = {ChannelMat.Channel.Name};
        ImaGIN_save_log(SpmFile, ['Convert: Selected channels     (' SpmFile ')'], ChanLabelsOut);
        % Save the list of channels in the output file
        ImaGIN_save_log(SpmFile, 'Convert: Removed channels:', sort(setdiff(ChanLabelsIn, ChanLabelsOut)));
    end
    set_final_status('OK');
end



%% ===== COPY FIELDS =====
% Copy fields from one structure to another.
%
% USAGE:  sDest = ImaGIN_copy_fields(sDest, sSrc, fieldNames) % Copy requested fields
%         sDest = ImaGIN_copy_fields(sDest, sSrc)             % Copy all fields
%
% INPUT:
%    - sDest : Destination structure 
%    - sSrc  : Source structure
%    - fieldNames : Cell array, list of fields to copy from sSrc to sDest
%
% Author: Francois Tadel, 2017
function sDest = ImaGIN_copy_fields(sDest, sSrc, fieldNames)
    if isempty(sSrc)
        return;
    end
    if (nargin < 3) || isempty(fieldNames)
        fieldNames = fieldnames(sSrc);
    end
    for i = 1:length(fieldNames)
        if isfield(sSrc, fieldNames{i})
            sDest.(fieldNames{i}) = sSrc.(fieldNames{i});
        end
    end
end
