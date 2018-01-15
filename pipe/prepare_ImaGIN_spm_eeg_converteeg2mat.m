function prepare_ImaGIN_spm_eeg_converteeg2mat(FileIn, FileOut)

% FileIn: path linking to the MEEG file to correct
% FileOut: output path


[Root,file,ext]=fileparts(FileIn);
clear S
S.FileOut=FileOut;
S.isSEEG=1;

switch lower(ext)
    
    case '.trc'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.channel=[];
        S.coarse=1;
        S.SaveFile=deblank(file);
        S.loadevents='yes';
        
    case '.msm'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.SEEG='Yes';
        S.coarse=1;
        S.SaveFile=deblank(file);
        
    case '.bin'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.SEEG='Yes';
        S.coarse=1;
        S.SaveFile=deblank(file);
        
    case '.asc'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';       
        S.Bipolar='No';
        S.coarse=1;
        S.SaveFile=deblank(file);
        S.Radc=PatientRadc;
        S.Nevent=PatientNevent;
        
    case '.edf'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.channel=[];
        S.SEEG='Yes';
        S.coarse=1;
        S.SizeMax=1e12;
        S.SaveFile=deblank(file);
        
        %     switch lower(tmp(end-1:end))
        %
    case '.e'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.channel=[];
        S.SEEG='Yes';
        S.coarse=1;
        S.SaveFile=['b' deblank(file)];
        
    case '.eeg'
        S.dataset=fullfile(FileIn);
        S.Atlas='Human';
        S.SEEG='Yes';
        S.coarse=1;
        S.SaveFile=['b' deblank(file)];
end

D = ImaGIN_spm_eeg_converteeg2mat(S);

end