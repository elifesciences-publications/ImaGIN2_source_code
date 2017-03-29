function sDest = ImaGIN_copy_fields(sDest, sSrc, fieldNames)
% Copy fields from one structure to another.
%
% USAGE:  sDest = ImaGIN_copy_fields(sDest, sSrc, fieldNames) % Copy requested fields
%         sDest = ImaGIN_copy_fields(sDest, sSrc)             % Copy all fields
%
% INPUT:
%    - sDest : Destination structure 
%    - sSrc  : Source structure
%    - fieldNames : Cell array, list of fields to copy from sSrc to sDest

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
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
%
% Author: Francois Tadel, 2017

if isempty(sSrc)
    return;
end

if (nargin < 3) || isempty(fieldNames)
    fieldNames = fieldnames(sSrc);
end

for i = 1:length(fieldNames)
    if isfield(sSrc, fieldNames{i})
        sDest.(fieldNames{i}) = sSrc.(fieldNames{i});
    end
end

