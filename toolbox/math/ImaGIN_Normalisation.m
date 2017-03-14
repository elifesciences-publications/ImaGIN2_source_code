function data=ImaGIN_Normalisation(data,dim,Baseline)

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

