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
% Copyright (c) 2000-2017 Inserm U1216
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
if isfield(S, 'SelectChannels') && ~isempty(S.SelectChannels)
    SelectChannels = S.SelectChannels;
else
    SelectChannels = [];
end

% Is output file defined
isOutputSet = isfield(S2, 'FileOut') && ~isempty(S2.FileOut);

% Processing starts
spm('Pointer','Watch');


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
    
    % Switch between file formats
    switch lower(fExt)
        case '.smr'
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'Bipolar', 'Montage', 'epochlength', 'coarse', 'channel'});
            D{i1} = ImaGIN_spm_eeg_rdata_spike2_mono(S2);
        
        case '.eeg'  % BrainAmp or Nihon Kohden
            % BrainAmp: There is a header in the same folder (.vhdr or .ahdr)
            if file_exist(fullfile(fPath, [fBase, '.vhdr'])) || file_exist(fullfile(fPath, [fBase, '.ahdr']))
                S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'Bipolar', 'Montage', 'epochlength', 'coarse', 'channel', 'SaveFile'});
                D{i1} = ImaGIN_spm_eeg_rdata_elan(S2);
            % Nihon Kohden
            else
                D{i1} = ImaGIN_convert_brainstorm(S2.Fdata, 'EEG-NK', [S2.FileOut, '.mat'], SelectChannels);
            end

        case {'.asc','.txt'}
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'coarse', 'Bipolar', 'channel', 'Radc', 'Nevent', 'filenamePos', 'filenameName', 'MontageName'});
            D{i1} = ImaGIN_spm_eeg_rdata_ascii(S2);

        case '.trc'
            % Old version: JPL & OD
            % S2 = ImaGIN_copy_fields(S2, S, {'pts', 'System', 'CreateTemplate', 'Fchannels', 'channel', 'coarse', 'Montage', 'MontageName', 'NeventType', 'event_file' , 'Bipolar', 'bipole', 'loadevents'});
            % D{i1} = ImaGIN_spm_eeg_rdata_micromed_mono(S2);

            % New version: Brainstorm
            D{i1} = ImaGIN_convert_brainstorm(S2.Fdata, 'EEG-MICROMED', [S2.FileOut, '.mat'], SelectChannels);

        case '.msm'
            S2 = ImaGIN_copy_fields(S2, S, {'Atlas', 'SEEG', 'Bipolar', 'coarse', 'SaveFile'});
            D{i1} = ImaGIN_spm_eeg_rdata_msm_mono(S2);

        case '.bin'
            S2 = ImaGIN_copy_fields(S2, S, {'Bipolar', 'Bipole', 'CreateTemplate', 'Montage', 'coarse', 'Nevent', 'SEEG', 'filenamePos', 'filenameName', 'MontageName', 'SaveFile'});
            D{i1} = ImaGIN_spm_eeg_rdata_deltamedbin_mono(S2);
       
        case '.edf'
            S2 = ImaGIN_copy_fields(S2, S, {'channel', 'Atlas', 'SEEG', 'Bipolar', 'coarse', 'SizeMax'});
            D = ImaGIN_spm_eeg_rdata_edf(S2);

        case '.e'
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Montage', 'coarse', 'SEEG', 'filenamePos', 'filenameName', 'MontageName', 'SaveFile', 'channel'});
            D = ImaGIN_spm_eeg_rdata_nicolet_mono(S2);

        otherwise
            error('Unknown format');
    end
end
% Processing stops
spm('Pointer','Arrow');
end



%% ===== CONVERT WITH BRAINSTORM I/O =====
% Converts EEG data from a native file to SPM-format, and saves a log of the channel selection
function D = ImaGIN_convert_brainstorm(InputFile, FileFormat, OutputFile, SelChannels) %#ok<STOUT>
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
        % NIHON KOHDEN .EEG
        case 'EEG-NK'
            [sFileIn, ChannelMat] = in_fopen_nk(InputFile);
        otherwise
            error(['Unsupported file format: ', FileFormat]);
    end
    % Read entire file with Brainstorm functions
    [F, TimeVector] = in_fread(sFileIn, ChannelMat, 1, []);
    
    % Save the original list of channels in a log file
    ChanLabelsIn = {ChannelMat.Channel.Name};
    ImaGIN_save_log(OutputFile, ['Available channels     (' InputFile ')'], ChanLabelsIn);

    % Auto-detect good SEEG channels
    if isempty(SelChannels)
        % Get channels classified as EEG
        iEEG = channel_find(ChannelMat.Channel, 'EEG,SEEG,ECOG');
        % Detect channels of interest
        iSelEeg = ImaGIN_select_channels({ChannelMat.Channel(iEEG).Name});
        % Convert indices back to the original list of channels
        iSel = iEEG(iSelEeg);
    % Selected channels are passed in input
    else
        % List of indices: use directly: [1 2 3 ...]
        if isnumeric(SelChannels)
            iSel = SelChannels;
        % Cell array of channel labels: {'A1','A2',...}
        elseif iscell(SelChannels)
            iSel = [];
            for i = 1:length(SelChannels)
                iSel = [iSel, find(strcmpi(lower(SelChannels{i}), lower(ChanLabelsIn)))];
            end
        % String: 'A1 A2 ...'
        elseif ischar(SelChannels)
            SelChannels = str_split(SelChannels, sprintf(' ,;\t'));
            iSel = [];
            for i = 1:length(SelChannels)
                iSel = [iSel, find(strcmpi(lower(SelChannels{i}), lower(ChanLabelsIn)))];
            end
        else
            error('Invalid value for parameter "SelectChannels".');
        end
    end
    % If no channels are selected
    if isempty(iSel)
        error('No valid channel names were found.');
    end
    % Keep only these ones in the data
    ChannelMat.Channel = ChannelMat.Channel(iSel);
    F = F(iSel,:);
    sFileIn.channelflag = sFileIn.channelflag(iSel);

    % Export to SPM format
    sFileOut = out_fopen_spm(OutputFile, sFileIn, ChannelMat);
    out_fwrite_spm(sFileOut, [], [], F);
    % Load new file to return the D structure
    load(OutputFile);

    % Save the list of channels in the output file
    ChanLabelsOut = {ChannelMat.Channel.Name};
    ImaGIN_save_log(OutputFile, ['Selected channels     (' OutputFile ')'], ChanLabelsOut);
    % Save the list of channels in the output file
    ImaGIN_save_log(OutputFile, 'Removed channels:', sort(setdiff(ChanLabelsIn, ChanLabelsOut)));
end



