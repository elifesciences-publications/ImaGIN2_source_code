function [H,D,N]=ImaGIN_permutation_entropy(data,k,TimeDelay,Prior)
%
% INPUTS:
%    - data: vector
%
% REFERENCES:
%    - Bandt & Pompe, 2002; Li X et al, 2007

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


Permutations=perms(1:k);
Nperm=size(Permutations,1);
N=zeros(length(TimeDelay),Nperm);
H=zeros(length(TimeDelay),1);
D=zeros(length(TimeDelay),1);
for i0=1:length(TimeDelay)
    NSample=length(data)-(k-1)*TimeDelay(i0);
    for i1=1:NSample
        index=i1:TimeDelay(i0):i1+(k-1)*TimeDelay(i0);
        [tmp1,tmp2]=sort(data(index));
        for i2=1:Nperm
            if isequal(tmp2,Permutations(i2,:))
                N(i0,i2)=N(i0,i2)+1;
                break
            end
        end
    end
    N(i0,:)=N(i0,:)./NSample;
    Index=find(N(i0,:)>0);
    H(i0)=-sum(N(i0,Index).*log(N(i0,Index)));    %permutation entropy
    H(i0)=log(factorial(k))./H(i0);             %normalisation factor (>1 if there exists a deterministic dynamics)
    if nargin==3
        D(i0)=sqrt(factorial(k)/(factorial(k)-1))*sqrt(sum((N(i0,:)-1/factorial(k)).^2));      %Dissimilarity with flat distribution (white noise)
    elseif nargin==4
        D(i0)=sqrt(factorial(k)/(factorial(k)-1))*sqrt(sum((N(i0,:)-transpose(Prior*pinv(Prior)*N(i0,:)')).^2));      %Dissimilarity with prior distribution of ordinal patterns (Prior)
    end
end
H=mean(H,1);
D=mean(D,1);
N=mean(N,1);




