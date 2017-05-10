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
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2017 Inserm U1216
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

% Remove all the channels with more than 18x the same tag (this is not SEEG)
if isSEEG
    uniqueTags = unique(chTags);
    for i = 1:length(uniqueTags)
        % Get channels of this tag
        iTag = find(strcmpi(uniqueTags{i}, chTags));
        % Remove if more than 18
        if (length(iTag) >= 18)
            chNames(iTag) = {'XXXXX'};
            chNamesClean(iTag) = {'XXXXX'};
            chTags(iTag) = {'XXXXX'};
        end
    end
end

% SEEG: Get the index of each channel
if isSEEG
    % Get indices
    chInd = cellfun(@(c)c(ismember(c, '0123456789')), chNamesClean, 'UniformOutput', 0);
    isNoInd = cellfun(@isempty, chInd);
    % Convert indices from string to integers
    chInd(~isNoInd) = cellfun(@str2num, chInd(~isNoInd), 'UniformOutput', 0);
    isNoInd = isNoInd | cellfun(@isempty, chInd);
% EEG: The separate name/index does not make sense
else
    chInd = num2cell(1:length(chNamesClean));
    isNoInd = zeros(1,length(chNamesClean));
end

% Process all the channels
iSel = [];
iEcg = [];
for i = 1:length(chNames)
    % No index or does not end with a digit
    if isSEEG && (isNoInd(i) || ~ismember(chNames{i}(end), '0123456789'))
        continue;
    % Does not start with a letter
    elseif isSEEG && ~ismember(lower(chNames{i}(1)), 'abcdefghijklmnopqrstuvwxyz')
        continue;
    % Unwanted labels
    elseif ismember(lower(chTags{i}), {'mark', 'dc', 'emg', 'eog', 'veo', 'heo', 'veog', 'heog', 'myo', 'oc', 'dd', 'dg', 'el', 'ref', 'eegref', 'eref', 'vref', 'ref', 'pulse', 'mast', 'spo2'})
        continue;
    % Unwanted labels
    elseif ~isempty(strfind(lower(chTags{i}), 'eog')) || ~isempty(strfind(lower(chTags{i}), 'ref'))
        continue;
    % Unwanted EEG labels
    elseif isSEEG && ismember(lower(chTags{i}), {'cz', 'fz', 'pz', 'oz', 'nz', 'fpz'})
        continue;
    % Unwanted FPx/Fx/Cx/Tx/Pz/Ox labels
    elseif isSEEG && ismember(lower(chNames{i}), {'fp1','fp2'}) && ~any(ismember({'fp3','fp4'}, lower(chNames)))
        continue;
    elseif isSEEG && ismember(lower(chNames{i}), {'f3','f4','f7','f8'}) && ~any(ismember({'f1','f2'}, lower(chNames)))
        continue;
    elseif isSEEG && ismember(lower(chNames{i}), {'c3','c4'}) && ~any(ismember({'c1','c2'}, lower(chNames)))
        continue;
    elseif isSEEG && ismember(lower(chNames{i}), {'t3','t4','t5','t6'}) && ~any(ismember({'t1','t2'}, lower(chNames)))
        continue;
    elseif isSEEG && ismember(lower(chNames{i}), {'p3','p4'}) && ~any(ismember({'p1','p2'}, lower(chNames)))
        continue;
    elseif isSEEG && ismember(lower(chNames{i}), {'o1','o2'}) && ~any(ismember({'o3','o4'}, lower(chNames)))
        continue;
    % ECG: Accept (should be labelled as such)
    elseif ismember(lower(chTags{i}), {'ecg', 'ekg'})
        iEcg(end+1) = i;
    % Otherwise: accept
    else
        iSel(end+1) = i;
    end
end

% Sort channels by tag and index
chTags = chTags(iSel);
chInd = chInd(iSel);
uniqueTags = unique(chTags);
iOrderSel = [];
for i = 1:length(uniqueTags)
    iTag = find(strcmpi(chTags, uniqueTags{i}));
    [tmp, iOrderInd] = sort([chInd{iTag}]);
    iOrderSel = [iOrderSel, iTag(iOrderInd)];
end

% Add ECG channels at the end of the file
iSel = [iSel(iOrderSel), iEcg];



