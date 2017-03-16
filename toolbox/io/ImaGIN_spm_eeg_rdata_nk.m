function D = ImaGIN_spm_eeg_rdata_nk(S)
% Converts EEG data from Nihon Kohden .EEG/.LOG/.PNT/.21E format to SPM-format
%
% USAGE: D = ImaGIN_spm_eeg_rdata_nk(S)
%
% Author: Francois Tadel, 2017

% Check fields in input
if (nargin == 0) || isempty(S) || ~isfield(S, 'Fdata') || ~isfield(S, 'FileOut')
    error(['Usage: D = ImaGIN_spm_eeg_rdata_nk(S)' 10 'Structure S must include fields: Fdata, FileOut']);
end

% Read Nihon Kohden file with Brainstorm functions
[sFileIn, ChannelMat] = in_fopen_nk(S.Fdata);
[F, TimeVector] = in_fread(sFileIn, ChannelMat, 1, []);

% TODO: FILTER CHANNELS OF INTEREST
warning('TODO: Filter channels of interest.')

% Output file
[fPath,fBase,fExt] = fileparts(S.Fdata);
OutputFile = fullfile(fPath, [fBase,'.mat']);
% Export to SPM format
sFileOut = out_fopen_spm(OutputFile, sFileIn, ChannelMat);
out_fwrite_spm(sFileOut, [], [], F);

% Load new file to return the D structure
load(OutputFile);







