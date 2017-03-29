function data=ImaGIN_Normalisation(data,dim,Baseline)
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
% Authors: Olivier David

if nargin==2||isempty(Baseline)
    Baseline=1:size(data,dim);
end
if dim==1
    mu=mean(data(Baseline,:),dim);
    sigma=std(data(Baseline,:),[],dim);
    data=(data-ones(size(data,1),1)*mu)./(ones(size(data,1),1)*sigma);
elseif dim==2
    mu=mean(data(:,Baseline),dim);
    sigma=std(data(:,Baseline),[],dim);
    data=(data-mu*ones(1,size(data,2)))./(sigma*ones(1,size(data,2)));
end

