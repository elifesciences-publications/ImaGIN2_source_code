function usSignal=ImaGIN_DownSample(Signal, fe, fs)
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

step1=round(fe/fs);
step2=ceil(fe/fs);

usSignal=zeros(step2,length(1:step1:length(Signal)));
for i=1:step2
%     s=Signal(1:step1:end);
%     i_sp1=floor(i/(fe/fs));
%     if length(s)>1
%         usSignal(i,:)=cat(1,s(i_sp1+1:end),s(1:i_sp1));
%     else
%         usSignal(i,:)=s;
%     end
%     Signal=cat(1,Signal(end), Signal(1:end-1));

    index=i+[0:step1:length(Signal)-i];
    usSignal(i,1:length(index))=Signal(index);
end

end