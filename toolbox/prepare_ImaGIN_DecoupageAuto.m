function prepare_ImaGIN_DecoupageAuto(NameStim,FileIn,FileOut)
 
% FileIn: path linking to the MEEG file to dissect

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
% Authors: ?

clear S;
S.dataset = fullfile(FileIn); 
S.DirFileOut = FileOut;

if ischar(NameStim)
    NameStim=str2num(NameStim);
end

if ~isempty(NameStim) && NameStim ~= 0
    S.StimName = NameStim;
else
    S.StimName = []; 
end
D = ImaGIN_DecoupageAuto(S);

