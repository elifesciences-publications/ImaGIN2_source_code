function prepare_ImaGIN_BadChannel(trainBase, FileIn, FileOut)
   
% FileIn: path linking to the MEEG file to extract bad channel indices
clear S;
S.dataset = fullfile(FileIn); 
S.FileOut = FileOut;
S.trainBase = trainBase; % directory where a trained model resides
fprintf('prepare_ImaGIN_BadChannel: S.dataset: %s \n', S.dataset)
fprintf('prepare_ImaGIN_BadChannel: S.FileOut: %s \n', S.FileOut)
fprintf('prepare_ImaGIN_BadChannel: S.trainBase: %s \n', S.trainBase)
D = ImaGIN_BadChannel(S);
