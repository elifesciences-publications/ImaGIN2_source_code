function CM=ImaGIN_ConnectivityMatrix(N)

ncouple	= N*(N-1)/2;
if ncouple>0
    CM	= zeros(2,ncouple);
    cou	= 0;
    for ii	= 1:N-1
        j	= ii+1;
        cou	= cou+1;
        CM(1,cou)	= ii;
        CM(2,cou)	= j;
        while j < N
            j	= j+1;
            cou	= cou+1;
            CM(1,cou)	= ii;
            CM(2,cou)	= j;
        end
    end
else
    CM=[];
end
return
