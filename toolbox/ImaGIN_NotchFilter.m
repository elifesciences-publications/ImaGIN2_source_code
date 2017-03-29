function ImaGIN_NotchFilter(P)
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
% Authors: Olivier David

try
    t=P.Fname;
catch
    F  = spm_figure('GetWin','Interactive');
    figure(F);clf
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

try
    nHz = P.Freq;
catch
    nHz = 50;
end

try
    FileOut=P.FileOut;
end


for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Data=D(:,:,:);
    for i1=1:D.nchannels
        for i2=1:D.ntrials
            Data(i1,:,i2)=ImaGIN_notch(squeeze(Data(i1,:,i2)),D.fsample,nHz,250);
        end
    end
    
    %Save as a newfile
    try
        D=clone(D,FileOut(i0), [D.nchannels D.nsamples D.ntrials]);
    catch
        D=clone(D,['n' D.fname], [D.nchannels D.nsamples D.ntrials]);
    end
    
    D(:,:,:)=Data;
    save(D);
end
