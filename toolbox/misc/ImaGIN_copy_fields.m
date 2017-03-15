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

