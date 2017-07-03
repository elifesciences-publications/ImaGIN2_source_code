function ImaGIN_LowPassFilter(P)
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
% Copyright (c) 2017 Inserm U1216
% =============================================================================-
%
% Authors: Viateur Tuyisenge & Olivier David

try
    t=P.LFname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.05,'DesignMethod','butter');

for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Data=D(:,:,:);
    for i1=1:D.nchannels
        for i2=1:D.ntrials
            Data(i1,:,i2)=filtfilt(d1, squeeze(Data(i1,:,i2)));
        end
    end

    %Save as a newfile
    D=clone(D,['lpf_' D.fname], [D.nchannels D.nsamples D.ntrials]);
    D(:,:,:)=Data;
    save(D);
end