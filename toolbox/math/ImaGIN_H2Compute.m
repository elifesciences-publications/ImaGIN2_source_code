function H2 = ImaGIN_H2Compute(x,y,ratio)
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

NTime=length(x);
if ratio>1
    Horizon=ratio;
else
    Horizon=round(NTime*ratio);
end
%H2(1,:): prediction of x from y
%H2(2,:): prediction of y from x
H2 = zeros(2,2*Horizon+1);
for i3 =-Horizon:Horizon
    tmp2 = ImaGIN_shift(y,i3);
    H2(1,i3+Horizon+1)=H2Compute(tmp2(Horizon+1:NTime-Horizon),x(Horizon+1:NTime-Horizon));
    H2(2,i3+Horizon+1)=H2Compute(x(Horizon+1:NTime-Horizon),tmp2(Horizon+1:NTime-Horizon));
%     tmp2=ImaGIN_shift(x,-i3);
end


function H2 = H2Compute(x,y)

%these mario chavez (p. 65)
%intervalles construits avec le meme nombre de points

x=reshape(x,prod(size(x)),1);
y=reshape(y,prod(size(x)),1);
x=x-min(x);
y=y-min(y);
Maxx = max(x);
Maxy = max(y);
N=length(x);
NBin = floor(log(N)/log(2));

[xsort,xorder]=sort(x);

p = zeros(1,NBin);
q = zeros(1,NBin);
for i1=1:NBin
%     if i1~=NBin
%         tmp = find(x<i1*Maxx/NBin&x>=(i1-1)*Maxx/NBin);
%     else
%         tmp = find(x<=i1*Maxx/NBin&x>=(i1-1)*Maxx/NBin);
%     end
%     q(i1) = mean(y(tmp));
%     p(i1) = Maxx/(2*NBin)+(i1-1)*Maxx/NBin;
    tmp = (i1-1)*floor(N/NBin)+[1:floor(N/NBin)];
    q(i1) = mean(y(xorder(tmp)));
    p(i1) = mean(xsort(tmp));
end
muxy = zeros(size(y));
for i1=1:length(muxy)
    for i2=1:NBin-1
        if x(i1)<=p(i2+1)
            tmp=i2;
            break
        elseif i2==NBin-1
            tmp=i2;
        end
    end
    muxy(i1)=(q(tmp+1)-q(tmp))*(x(i1)-p(tmp))/(p(tmp+1)-p(tmp))+q(tmp);
end
H2 = 1-var(y-muxy)/var(y);
