function X=ImaGIN_FT1Surrogate(x,NRealisation)
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

N=length(x);

x=reshape(x,1,N);
X=zeros(NRealisation,N);

fx=fft(x);
if ceil(N/2)~=N/2
    N1=ceil(N/2);
    N2=ceil(N/2)-1;
else
    N1=ceil(N/2)+1;
    N2=ceil(N/2);
end
if mod(N,2)==1
    tmp1=2*pi*rand(NRealisation,ceil(N/2));
    tmp1=[tmp1 -tmp1(:,[ceil(N/2):-1:2])];
else
    tmp1=2*pi*rand(NRealisation,ceil(N/2)+1);
    tmp1=[tmp1 -tmp1(:,[ceil(N/2):-1:2])];
end
for i1=1:NRealisation
    X(i1,:)=real(ifft(abs(fx).*exp(1i.*(angle(fx)+tmp1(i1,:)))));
end

