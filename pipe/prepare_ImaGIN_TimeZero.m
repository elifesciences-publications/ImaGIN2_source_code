function prepare_ImaGIN_TimeZero(EventRef, Offset, FileIn,  FileOut)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% EventRef: type of event (default:1)
% Offset: set the time 0 at 0+Offset (default:0)

S.Filename=FileIn;

S.EventRef=EventRef;

if ischar(Offset)
    Offset=str2num(Offset);
end
S.Offset=Offset;

S.FileOut=FileOut;
ImaGIN_TimeZero(S);

end