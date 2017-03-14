function prepare_ImaGIN_BadChannel(FileIn,trainBase,FileOut)

% FileIn: path linking to the MEEG file to extract bad channel indices
clear S;
S.dataset = fullfile(FileIn); 
S.DirFileOut = FileOut;
S.trainBase = trainBase; % directory where a trained model resides
D = ImaGIN_BadChannel(S);