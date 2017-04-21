function D = ImaGIN_convert_nk(InputFile, OutputFile)
% Converts EEG data from Nihon Kohden .EEG/.LOG/.PNT/.21E format to SPM-format
%
% USAGE: D = ImaGIN_convert_nk(InputFile, OutputFile)

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
% Author: Francois Tadel, 2017


% Read Nihon Kohden file with Brainstorm functions
[sFileIn, ChannelMat] = in_fopen_nk(InputFile);
[F, TimeVector] = in_fread(sFileIn, ChannelMat, 1, []);

% Get channels classified as EEG
iEEG = channel_find(ChannelMat.Channel, 'EEG,SEEG,ECOG');
% Detect channels of interest
iSel = ImaGIN_select_channels({ChannelMat.Channel(iEEG).Name});
if isempty(iSel)
    error('No valid channel names were found.');
end
% Keep only these ones in the data
ChannelMat.Channel = ChannelMat.Channel(iEEG(iSel));
F = F(iEEG(iSel),:);
sFileIn.channelflag = sFileIn.channelflag(iEEG(iSel));

% Export to SPM format
sFileOut = out_fopen_spm(OutputFile, sFileIn, ChannelMat);
out_fwrite_spm(sFileOut, [], [], F);

% Load new file to return the D structure
load(OutputFile);







