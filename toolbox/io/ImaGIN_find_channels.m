function iChannels = ImaGIN_find_channels(ChannelLabels, SelChannels)
% Find channel indices from a list.

% INPUTS:
%    - ChannelLabels : Cell array of strings, list of available channel names.
%    - SelChannels   : Channels to select from this list, possible options:
%         |            - Array of indices ([1,2,...]), 
%         |            - Array of cells ({'A1','A2',...}) 
%         |            - String, channel names separated with comas or spaces ('A1 A2 ...')
%         |            - String, channel indices separated with spaces ('1 2 ...')

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

% List of indices: use directly: [1 2 3 ...]
if isnumeric(SelChannels)
    iChannels = SelChannels;
% Cell array of channel labels: {'A1','A2',...}
elseif iscell(SelChannels)
    iChannels = [];
    for i = 1:length(SelChannels)
        iChannels = [iChannels, find(strcmpi(lower(SelChannels{i}), lower(ChannelLabels)))];
    end
% String: 'A1 A2 ...' or '1 2 ...'
elseif ischar(SelChannels)
    % Try to read as channel names
    splitChannels = str_split(SelChannels, sprintf(' ,;\t'));
    iChannels = [];
    for i = 1:length(splitChannels)
        iChannels = [iChannels, find(strcmpi(lower(splitChannels{i}), lower(ChannelLabels)))];
    end
    % If nothing is found, try to read as channel indices
    if isempty(iChannels)
        iChannels = str2num(SelChannels);
    end
    % If nothing is found, try to interpret as a Matlab expression
    if isempty(iChannels)
        try
            iChannels = eval(SelChannels);
        catch
            iChannels = [];
        end
    end
else
    error('Invalid value for parameter SelChannels.');
end

% Remove channels that are not in the available range
iChannels(iChannels > length(ChannelLabels)) = [];
    
    


