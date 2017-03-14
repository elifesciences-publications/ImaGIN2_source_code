function AC = ImaGIN_Crosscorrelation(x,y,ratio)

x=x-mean(x);
y=y-mean(y);

NTime=length(x);
Horizon=round(NTime*ratio);
AC = zeros(1,2*Horizon+1);
for i3 =-Horizon:Horizon
    tmp2 = ImaGIN_shift(y,i3);
    S=sqrt(sum(x(Horizon+1:NTime-Horizon).^2).*sum(tmp2(Horizon+1:NTime-Horizon).^2));
    AC(i3+Horizon+1)=sum(x(Horizon+1:NTime-Horizon).*tmp2(Horizon+1:NTime-Horizon))./S;
end


