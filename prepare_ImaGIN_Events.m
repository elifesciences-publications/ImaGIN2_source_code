function prepare_ImaGIN_Events(Action, Nevent, EventName, FileDataIn, FileEventIn, FileOut)

% FileIn: full path linking to the MEEG file to correct
% DirStimIn: full path linking to the output directory of ImaGIN_StimDetect
% ('_StimulationIndex.txt' & '_Stimulation')
% DirOut: output path
% Action: 'Add' or 'Remove' an event
% Nevent: number of the event (to be registered in the MEEG file)
% (default:1)
% EventName: name of the event (default:'Stim')

S.Fname=FileDataIn;

S.Action=Action;
if ischar(Nevent)
    Nevent=str2num(Nevent);
end

S.Nevent=Nevent;
S.EventName{1}=EventName;
S.EventFileName{1}=FileEventIn;
S.FileOut=FileOut;
ImaGIN_Events(S);

close all

end