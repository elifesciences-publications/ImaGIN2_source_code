function prepare_ImaGIN_spm_eeg_converteeg2mat(FileIn, FileOut)

% FileIn: path linking to the MEEG file to correct
% FileOut: output path

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2017 Inserm U1216
% =============================================================================-
%
% Authors: ?

[Root,file,ext]=fileparts(FileIn);
clear S
S.FileOut=FileOut;

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