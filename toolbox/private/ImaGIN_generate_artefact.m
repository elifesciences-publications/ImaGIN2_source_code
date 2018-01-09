function y = ImaGIN_generate_artefact(tau,fe, art_duration, break_duration, total_window, mode)
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

i=zeros(1,round(total_window*fe));
i(1)=0;
if strcmp(mode, 'biphasic')
    i(2:round(art_duration*fe/2)+1)=1;
    i(round(art_duration*fe/2)+2:round((art_duration+break_duration)*fe/2)+2)=0;
    i(round((art_duration+break_duration)*fe/2)+3:round((art_duration+break_duration)*fe)+3)=-1;
elseif strcmp(mode, 'monophasic')
    i(2:round(art_duration*fe))=1;
end


MY=zeros(length(i),length(tau));
dt=1/fe;

for jj=1:length(i)
    if jj==1
        y = 0;
    end
    dy = (-y+i(jj))/tau;
    y = y + dy*dt;
    MY(jj,:)=y;
end

y=(i'-MY);
end






















