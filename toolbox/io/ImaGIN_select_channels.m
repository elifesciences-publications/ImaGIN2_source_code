function [iSel, iEcg] = ImaGIN_select_channels(chNames, isSEEG)
% IMAGIN_SELECT_CHANNELS Keep only channels of interest.
%
% USAGE:  [iSel, iEcg] = ImaGIN_select_channels(chNames, isSEEG=1)
%
% INPUT: 
%    - chNames : Cell-array of strings
%    - isSEEG  : 1 if the data is SEEG, 0 if regular EEG
%
% OUTPUT:
%    - iSel : Array of indices of the channels that are considered as valid
%    - iEcg : Array of indices of the channels that are identified as ECG

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% Copyright (c)2000-2017 Inserm
% =============================================================================-
%
% Authors: Francois Tadel, 2017

% Parse inputs
if (nargin < 2) || isempty(isSEEG)
    isSEEG = 1;
end

% Get all names: remove special characters
chNamesClean = cellfun(@(c)c(~ismember(c, ' .,?!-_@#$%^&*+*=()[]{}|/')), chNames, 'UniformOutput', 0);
% Separate characters and numbers in the names
chTags = cellfun(@(c)c(~ismember(c, '0123456789')), chNamesClean, 'UniformOutput', 0);

% Remove all the channels with more than 18x the same tag (not SEEG)
if isSEEG
    uniqueTags = unique(chTags);
    for i = 1:length(uniqueTags)
        % Get channels of this tag
        iTag = find(strcmpi(uniqueTags{i}, chTags));
        % Remove if more than 18
        if (length(iTag) >= 18)
            chNames(iTag) = [];
            chNamesClean(iTag) = [];
            chTags(iTag) = [];
        end
    end
end

% Get indices
chInd = cellfun(@(c)c(ismember(c, '0123456789')), chNamesClean, 'UniformOutput', 0);
isNoInd = cellfun(@isempty, chInd);
% chInd(~isNoInd) = cellfun(@str2num, chInd(~isNoInd));


% Process all the channels
iSel = [];
iEcg = [];
for i = 1:length(chNames)
    % No index or does not end with a digit
    if isNoInd(i) || ~ismember(chNames{i}(end), '0123456789')
        continue;
    % Does not start with a letter
    elseif ~ismember(lower(chNames{i}(1)), 'abcdefghijklmopqrstuvwxyz')
        continue;
    % Unwanted labels
    elseif ismember(lower(chTags{i}), {'mark', 'dc', 'emg', 'eog', 'veo', 'heo', 'veog', 'heog', 'myo', 'oc', 'dd', 'dg', 'el', 'ref', 'eegref', 'eref', 'vref', 'ref', 'pulse', 'mast'})
        continue;
    % Unwanted EEG labels
    elseif isSEEG && ismember(lower(chTags{i}), {'cz', 'fz', 'pz', 'oz', 'nz', 'fpz'})
        continue;
    % ECG: Accept (should be labelled as such)
    elseif ismember(lower(chTags{i}), {'ecg', 'ekg'})
        iSel(end+1) = i;
        iEcg(end+1) = i;
    % Otherwise: accept
    else
        iSel(end+1) = i;
    end
end


