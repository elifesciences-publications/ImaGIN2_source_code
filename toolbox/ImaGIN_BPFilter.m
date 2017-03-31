function ImaGIN_BPFilter(P)
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
    t = P.Fname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

F  = spm_figure('GetWin','Interactive');
figure(F);clf
try
    BP = P.BP;
catch
    BP = spm_input('Bandpass filter (e.g. 5 40)','1','r');
end

for i0=1:size(t,1)
    T = deblank(t(i0,:));
    D = spm_eeg_load(T);
    Data = D(:,:,:);
    for i1 = 1:D.nchannels
        for i2 = 1:D.ntrials 
            Data(i1,:,i2) = ImaGIN_bandpassFilter(squeeze(Data(i1,:,i2)),D.fsample,BP(1),BP(2));
        end
    end

    % Save as a newfile
    Prefix = [num2str(BP(1)) '-' num2str(BP(2))];
    Dnew = clone(D, [Prefix '_' fname(D)], [D.nchannels D.nsamples D.ntrials]);
    Dnew(:,:,:) = Data;
    save(Dnew);
end
