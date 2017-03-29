function prepare_bad_channels(Badp, FileIn, FileOut)
% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% Badp: path linking to the file with the name of bad channels

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
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
%
% Authors: ?

Bad=load(Badp);
[DirOut,NameOut,~]=fileparts(FileOut);
if ~isempty(Bad)
    S.D=FileIn;
    D=spm_eeg_load(S.D);
    D=badchannels(D, Bad, 1);
    nameDproc=fullfile(DirOut, NameOut);
    D2=clone(D,nameDproc, [D.nchannels D.nsamples D.ntrials]);
    D2(:,:,:)=D(:,:,:);
    save(D);
end

end