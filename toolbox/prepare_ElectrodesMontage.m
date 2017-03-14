function prepare_ElectrodesMontage(PatientName, RootImplant, FileIn, DirOut)

[Root,file,ext]=fileparts(FileIn);
% tmp=deblank([file '.' ext]);

switch lower(ext)
    
    case '.e'
        %Add electrode position
        clear S
        S.Fname=fullfile(Root,PatientName,[spm_str_manip(file,'r') '.mat']);
        S.filenamePos=fullfile(RootImplant,PatientName,'Implantation','Electrodes_Pos_MNI_Bip.txt');
        S.filenameName=fullfile(RootImplant,PatientName,'Implantation','Electrodes_Name_Bip.txt');
        S.Out=DirOut;
        D = ImaGIN_Electrode(S);
        S.SaveFile=deblank(file);
        disp(S.SaveFile(1:end-2));
        
    otherwise
        %Add electrode position
        clear S
        S.Fname=fullfile(Root,PatientName,[spm_str_manip(file,'r') '.mat']);
        S.filenamePos=fullfile(RootImplant,PatientName,'Implantation','Electrodes_Pos_MNI.txt');
        S.filenameName=fullfile(RootImplant,PatientName,'Implantation','Electrodes_Name.txt');
        S.Out=DirOut;
        D = ImaGIN_Electrode(S);
        
        %Longitudinal bipolar montage
        clear S
        S.Fname=fullfile(Root,PatientName,[spm_str_manip(file,'r') '.mat']);
        S.SaveFile=['b' spm_str_manip(file,'r')];
        S.Out=DirOut;
        D = ImaGIN_BipolarMontage(S);
end

S.SaveFile=deblank(file);
disp(S.SaveFile(1:end-length(ext)-1));

end