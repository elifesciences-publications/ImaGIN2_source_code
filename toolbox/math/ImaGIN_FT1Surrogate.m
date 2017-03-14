function X=ImaGIN_FT1Surrogate(x,NRealisation)

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

