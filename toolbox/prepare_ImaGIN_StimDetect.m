function prepare_ImaGIN_StimDetect(StimContinuous, FileIn, FileOut)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% Stimstart: beginning of the stimulation (time in s - default:[])
% StimEnd: end of the stimulation (time in s - default:[])
% StimContinuous: indiquates whether the stimulation is continuous (1) or
% not (0)

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

S.Filename=FileIn;
S.StimFreq=1;

if strcmp(StimContinuous,'False')
    StimContinuous=false;
else if strcmp(StimContinuous,'True')
        StimContinuous=true;
    end
end

try
    S.StimContinuous=StimContinuous;
catch
    S.StimContinuous=0;
end

S.FileOut=FileOut;

S.Channels=[];
ImaGIN_StimDetect(S);


end