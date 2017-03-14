
function prepare_BipolarMontage(Save,FileIn, FileOut)

[Root,file,ext]=fileparts(FileIn);
%Longitudinal bipolar montage
clear S
S.Filename=FileIn;
S.FileOut=FileOut;

if strcmp(Save,'False')
    Save=false;
else if strcmp(Save,'True')
        Save=true;
    end
end

S.Save=Save;
D = ImaGIN_BipolarMontage(S);

S.SaveFile=deblank(file);
disp(S.SaveFile(1:end-length(ext)-1));

end