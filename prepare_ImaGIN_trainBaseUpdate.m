function prepare_ImaGIN_trainBaseUpdate(inDir)
% FileIn: path linking to the Cropped MEEG object 
clear S;
S.DirIn = inDir;
ImaGIN_trainBaseUpdate(S)
