function D = ImaGIN_convert_micromed(InputFile, OutputFile)
% Converts EEG data from a Micromed .TRC file to SPM-format
%
% USAGE: D = ImaGIN_convert_micromed(InputFile, OutputFile)

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
[sFileIn, ChannelMat] = in_fopen_micromed(InputFile);
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

% Export to SPM format
sFileOut = out_fopen_spm(OutputFile, sFileIn, ChannelMat);
out_fwrite_spm(sFileOut, [], [], F);

% Load new file to return the D structure
load(OutputFile);







