function prepare_ImaGIN_Validate_StimNames(defaultPulseDuration, patientName, FileIn, FileOut)
clear S;
S.dataset = fullfile(FileIn); 
S.DirFileOut = FileOut;
S.defaultPulseDuration = defaultPulseDuration;
S.patientName = patientName;
ImaGIN_Validate_StimNames(S);

close all

end
