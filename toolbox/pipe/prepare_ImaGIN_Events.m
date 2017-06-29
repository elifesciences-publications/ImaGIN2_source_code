function prepare_ImaGIN_Events(Action, Nevent, EventName, FileDataIn, FileEventIn, FileOut)

% FileIn: full path linking to the MEEG file to correct
% DirStimIn: full path linking to the output directory of ImaGIN_StimDetect
% ('_StimulationIndex.txt' & '_Stimulation')
% DirOut: output path
% Action: 'Add' or 'Remove' an event
% Nevent: number of the event (to be registered in the MEEG file)
% (default:1)
% EventName: name of the event (default:'Stim')

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

S.Fname = FileDataIn;

S.Action=Action;
if ischar(Nevent)
    Nevent=str2num(Nevent);
end

S.Nevent=Nevent;
S.EventName{1}=EventName;
S.EventFileName{1}=FileEventIn;
S.FileOut=FileOut;
ImaGIN_Events(S);

end