function prepare_ImaGIN_CropAuto(NameStim,FileIn,FileOut)
   
% FileIn: path linking to the MEEG file to dissect

clear S;
S.dataset = fullfile(FileIn); 
S.DirFileOut = FileOut;

if ischar(NameStim)
    NameStim=str2num(NameStim);
end

if ~isempty(NameStim) && NameStim ~= 0
    S.StimName = NameStim;
else
    S.StimName = []; 
end
D = ImaGIN_CropAuto(S);

close all 

end
