function prepare_ImaGIN_Validate_StimNames(defaultPulseDuration, FileIn, FileOut)
clear S;
S.dataset = fullfile(FileIn); 
S.DirFileOut = FileOut;
S.defaultPulseDuration = defaultPulseDuration;
ImaGIN_Validate_StimNames(S);