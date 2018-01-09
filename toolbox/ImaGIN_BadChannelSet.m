function D = ImaGIN_BadChannelSet(S)
% Set the list of bad channels for a SPM .mat/.dat file.

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
% Authors: Francois Tadel, 2017


%% ===== GET INPUTS =====
% Nothing in input
if (nargin < 1) || isempty(S)
    S = struct();
end
% Get filename
if isfield(S, 'Fname') && ~isempty(S.Fname)
    Filename = S.Fname;
else
    Filename = spm_select(1, '\.mat$', 'Select data file');
end
if isempty(Filename)
    D = [];
    return;
end
% Get bad channels
if isfield(S, 'BadChannels') && ~isempty(S.BadChannels)
    BadChannels = S.BadChannels;
else
    BadChannels = spm_input('Bad channels (names or indices)', '+1', 's');
end


%% ===== UPDATE FILE =====
% Load input file
D = spm_eeg_load(Filename);
% Get channel indices
ChanLabels = chanlabels(D);
iBadChan = ImaGIN_find_channels(ChanLabels, BadChannels);
% No bad channels: nothing else to do
if isempty(iBadChan)
    disp(['ImaGIN> Warning: No channel found for input "' BadChannels '"']);
    return;
end
% Set bad channels
D = badchannels(D, iBadChan, 1);
% Update file
save(D);

% Add the list of bad channels to the logs
ImaGIN_save_log(fullfile(D), 'Bad channels: ', ChanLabels(iBadChan));


