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
AllNames = cellfun(@(c)c(~ismember(c, ' .,?!-_@#$%^&*+*=()[]{}|/')), chNames, 'UniformOutput', 0);
AllTags  = cell(size(AllNames));
AllInd   = cell(size(AllNames));
isNoInd  = zeros(size(AllNames));
% Separate characters and numbers in the names
for i = 1:length(AllNames)
    % Find the last letter in the name
    iLastLetter = find(~ismember(AllNames{i}, '0123456789'), 1, 'last');
    if isempty(iLastLetter)
        iLastLetter = 0;
    end
    AllTags{i} = AllNames{i}(1:iLastLetter);
    % For EEG: The separate name/index does not make sense
    if ~isSEEG
        AllInd{i} = sprintf('%d',i);
        isNoInd(i) = 0;
    % If there are digits at the end of the name: use them as the index of the contact
    elseif (iLastLetter < length(AllNames{i}))
        AllInd{i} = AllNames{i}(iLastLetter+1:end);
    else
        isNoInd(i) = 1;
        AllInd{i} = '0';
    end
end
% Convert indices to double values
AllInd = cellfun(@str2num, AllInd);
    
% Remove all the channels with more than 18x the same tag (this is not SEEG)
if isSEEG
    uniqueTags = unique(AllTags);
    for i = 1:length(uniqueTags)
        % Get channels of this tag
        iTag = find(strcmpi(uniqueTags{i}, AllTags));
        % Remove if more than 18 (except for Salpetriere electrodes with digits in the name)
        if ((length(iTag) > 18) && (any(iTag < 10) || any(iTag >= 30)))
            AllNames(iTag) = {'XXXXX'};
            AllTags(iTag) = {'XXXXX'};
        end
        % Removed: ((length(iTag) < 2) && ~any(chNames{iTag(1)} == '-'))   % less than 2 (but only if not designating a bipolar montage, with '-' in the name) 
        % Keep electrodes with one contact only
    end
end

% Process all the channels
iSel = [];
iEcg = [];
for i = 1:length(AllNames)
    % No index or does not end with a digit
    if isSEEG && (isNoInd(i) || ~ismember(AllNames{i}(end), '0123456789'))
        continue;
    % Does not contain at least a letter
    elseif isSEEG && ~any(ismember(lower(AllNames{i}), 'abcdefghijklmnopqrstuvwxyz'))
        continue;
    % Unwanted labels
    elseif ismember(lower(AllTags{i}), {'xxxxx', 'mark', 'dc', 'emg', 'eog', 'veo', 'heo', 'veog', 'heog', 'myo', 'oc', 'dd', 'dg', 'el', 'ref', 'eegref', 'eref', 'vref', 'ref', 'pulse', 'mast', 'spo2'})
        continue;
    % Unwanted labels
    elseif ~isempty(strfind(lower(AllTags{i}), 'eog')) || ~isempty(strfind(lower(AllTags{i}), 'ref'))
        continue;
    % Unwanted EEG labels
    elseif isSEEG && ismember(lower(AllTags{i}), {'cz', 'fz', 'pz', 'oz', 'nz', 'fpz'})
        continue;
    % Unwanted FPx/Fx/Cx/Tx/Pz/Ox labels
    elseif isSEEG && ismember(lower(AllNames{i}), {'fp1','fp2'}) && ~any(ismember({'fp3','fp4'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'f3','f4','f7','f8'}) && ~any(ismember({'f1','f2'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'ft1','ft2','ft9','ft10'}) && ~any(ismember({'ft5','ft6'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'tp9','tp10'}) && ~any(ismember({'tp1','tp2'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'c3','c4'}) && ~any(ismember({'c1','c2'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'t3','t4','t5','t6','t9','t10'}) && ~any(ismember({'t1','t2'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'p3','p4'}) && ~any(ismember({'p1','p2'}, lower(AllNames)))
        continue;
    elseif isSEEG && ismember(lower(AllNames{i}), {'o1','o2'}) && ~any(ismember({'o3','o4'}, lower(AllNames)))
        continue;
    % ECG: Accept (should be labelled as such)
    elseif ismember(lower(AllTags{i}), {'ecg', 'ekg'})
        iEcg(end+1) = i;
    % Otherwise: accept
    else
        iSel(end+1) = i;
    end
end

% Sort channels by tag and index
AllTags = AllTags(iSel);
AllInd = AllInd(iSel);
uniqueTags = unique(AllTags);
iOrderSel = [];
for i = 1:length(uniqueTags)
    iTag = find(strcmpi(AllTags, uniqueTags{i}));
    [tmp, iOrderInd] = sort(AllInd(iTag));
    iOrderSel = [iOrderSel, iTag(iOrderInd)];
end

% Add ECG channels at the end of the file
iSel = [iSel(iOrderSel), iEcg];



