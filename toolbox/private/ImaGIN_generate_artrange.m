function mArt = ImaGIN_generate_artrange(rcmin, rcmax, nb, fe, art_duration, break_duration, window, mode)
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
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David

id='MATLAB:NonIntegerInput';
warning('off',id);

prange=linspace(rcmin,rcmax,nb);

mArt=zeros(round(window*fe),nb);
for ii=1:length(prange)
    tau=prange(ii);
    mArt(:,ii)=ImaGIN_generate_artefact(tau,fe,art_duration, break_duration, window, mode);
end

end


