function prepare_ImaGIN_StimDetect(StimContinuous, FileIn, FileOut)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% Stimstart: beginning of the stimulation (time in s - default:[])
% StimEnd: end of the stimulation (time in s - default:[])
% StimContinuous: indiquates whether the stimulation is continuous (1) or
% not (0)

S.Fname=FileIn;
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