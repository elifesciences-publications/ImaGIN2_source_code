function D = ImaGIN_spm_eeg_converteeg2mat(S)
% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_converteeg2mat(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be either 'CNT' or 'BDF'
%    Mname     - char matrix of input file name(s)
%    Fchannels - String containing name of channel template file
%_______________________________________________________________________
%
% spm_eeg_converteeg2mat is a user interface to convert EEG-files from their
% native format to SPM's data format. This function assembles some
% necessary information before branching to the format-specific conversion
% routines.
% The user has to specify, by either using struct S or the GUI, a 'channel
% template file' that contains information about the (approximate) spatial
% positions of the channels.
%
% Output: The converted data are written to files. The header
% structs, but not the data, are returned in D as a cell vector of structs.
%_______________________________________________________________________
%
% Additional formats can be added by (i) extending the code below in a
% straightforward fashion, (ii) providing a new channel template file and
% (iii) adding the actual conversion routine to the SPM-code.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_converteeg2mat.m 193 2005-07-01 13:37:59Z james $

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

% Ask species (if not defined in input)
if ~isfield(S2, 'Atlas')
    S2.Atlas = spm_input('Select atlas', '+1','Human|Rat|Mouse');
end

% Get output file (if not defined in input)
if ~isfield(S2, 'FileOut')
    S2.FileOut = spm_str_manip(S.dataset(1:end-4),'t');
end

% Processing starts
spm('Pointer','Watch');


% Loop on input files
Nfiles = size(S.dataset, 1);
D = cell(1,Nfiles);
for i1 = 1:Nfiles
    % Get file extension
    S2.Fdata = deblank(S.dataset(i1, :));
    [fPath, fBase, fExt] = fileparts(S2.Fdata);
               
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
                D{i1} = ImaGIN_spm_eeg_rdata_nk(S2);
            end

        case {'.asc','.txt'}
            S2 = ImaGIN_copy_fields(S2, S, {'CreateTemplate', 'Fchannels', 'coarse', 'Bipolar', 'channel', 'Radc', 'Nevent', 'filenamePos', 'filenameName', 'MontageName'});
            D{i1} = ImaGIN_spm_eeg_rdata_ascii(S2);

        case '.trc'
            S2 = ImaGIN_copy_fields(S2, S, {'pts', 'System', 'CreateTemplate', 'Fchannels', 'channel', 'coarse', 'Montage', 'MontageName', 'NeventType', 'event_file' , 'Bipolar', 'bipole', 'loadevents'});
            D{i1} = ImaGIN_spm_eeg_rdata_micromed_mono(S2);

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




