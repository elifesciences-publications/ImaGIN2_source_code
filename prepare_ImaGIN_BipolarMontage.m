function prepare_ImaGIN_BipolarMontage(Save,FileIn, FileOut)
   
[Root,file,ext]=fileparts(FileIn);
%Longitudinal bipolar montage
clear S
S.Fname=FileIn;
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

set_final_status('OK')

close all

end
