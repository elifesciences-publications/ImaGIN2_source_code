function prepare_ImaGIN_SpikesDetection(FileIn, FileOut)
   
% FileIn: path linking to the SPM object 
clear S;
S.dataset = fullfile(FileIn); 
S.FileOut = FileOut;
fprintf('prepare_ImaGIN_SpikesDetection: S.dataset: %s \n', S.dataset)
fprintf('prepare_ImaGIN_SpikesDetection: S.FileOut: %s \n', S.FileOut)
D = ImaGIN_SpikesDetection(S);
end