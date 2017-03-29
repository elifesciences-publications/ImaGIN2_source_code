function prepare_ImaGIN_TimeZero(EventRef, Offset, FileIn,  FileOut)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% EventRef: type of event (default:1)
% Offset: set the time 0 at 0+Offset (default:0)

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
% Authors: ?

S.Filename=FileIn;

S.EventRef=EventRef;

if ischar(Offset)
    Offset=str2num(Offset);
end
S.Offset=Offset;

S.FileOut=FileOut;
ImaGIN_TimeZero(S);

end