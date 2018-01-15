function prepare_ImaGIN_BadChannelQualityControl(FileIn)
% FileIn: path linking to the Cropped MEEG object 
clear S;
S.dataset = fullfile(FileIn); 
D = ImaGIN_BadChannelQualityControl(S);
