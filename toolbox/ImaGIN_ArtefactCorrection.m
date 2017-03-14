function D=ImaGIN_ArtefactCorrection(S)

% Correct artefacts in frequency space
% Ref: Ramirez & Baillet, 2011

try
    t=S.Filename;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end
for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Time=time(D);
    
    F  = spm_figure('GetWin','Interactive');
    figure(F);clf
    
    
    try
        EventType = S.EventType;
    catch
        EventType=spm_input('Events of interest', '+1', 'i',1);
    end
    try
        StartInterpolation = S.StartInterpolation;
    catch
        StartInterpolation=spm_input('Start of interpolation [sec]', '+1', 'r',-0.012);
    end
    try
        EndInterpolation = S.EndInterpolation;
    catch
        EndInterpolation=spm_input('End of interpolation [sec]', '+1', 'r',0.012);
    end
    
    int_corr=int2str(EndInterpolation*10^3);
    StartInterpolation=round(StartInterpolation*fsample(D));
    EndInterpolation=round(EndInterpolation*fsample(D));
    
    Bad=badchannels(D);
    Nchan=nchannels(D);
    Good=setdiff(1:Nchan,Bad);
    Nchan=length(Good);
    
    ev=events(D);
    
    Data=D(:,:);  
    
    NCompo=2;
    ind12=[];
    ind22=[];
    for i1=1:length(ev)
        if ~isempty(intersect(EventType,ev(i1).value))
            ind=indsample(D,ev(i1).time);
            ind12=[ind12 ind+StartInterpolation];
            ind22=[ind22 ind+EndInterpolation];
        end
    end
    
    ind12=sort(ind12);
    ind22=sort(ind22);
    Index2=[];
    
    for i1=1:length(ind12)
        Index2=[Index2 ind12(i1):ind22(i1)];
    end
    
    Index2=Index2(~isnan(Index2));
    
    Datatmp=Data(:,Index2);
    M=zeros(size(Datatmp));
    
    for i1=1:length(ind12)
        M(:,(i1-1)*length(Index2)/length(ind12)+[1:length(Index2)/length(ind12)])=...
            mean(Datatmp(:,(i1-1)*length(Index2)/length(ind12)+[1:length(Index2)/length(ind12)]),2)*ones(1,length(Index2)/length(ind12));
    end
    
    data=transpose(Data(:,Index2)-M);
    [u s v]=svd(data,0);
    d2=[zeros(size(data,1),1) abs(diff(data,2,2)) zeros(size(data,1),1)];
    d2=mean(d2,2);
    d2=ImaGIN_Normalisation(d2',2);
    d2=d2+abs(min(d2));
    sizeWindow=EndInterpolation-StartInterpolation+1; 

    
if D.fsample==512
 model=convn(abs(u(:,1)),ones(1,2)');
elseif D.fsample==1024
    model=convn(abs(u(:,1)),ones(1,3)');
else model=abs(u(:,1));
end

 c=FindClusters(model,sizeWindow);
    Index=c(find(c));
    Index=sort(Index);
    
    f=abs(fft(u));
    CC=zeros(1,size(f,2));
    
    for i1=1:size(u,2)
        cc=corrcoef(f(Index,1),f(Index,i1));
        CC(1,i1)=abs(cc(2));
    end
    
    for i1=1:size(u,2)
        cc=corrcoef(u(Index,1),u(Index,i1));
        CC(2,i1)=abs(cc(2));
    end
   
   CC=mean(CC,1);
   
   if graythresh(CC)<0.5
  CCthresh=graythresh(CC);
   else CCthresh=0.15;
   end
   
    CompoArtefact=find(abs(CC)>CCthresh);
    Time3=Time(Index2);
    Index3=setdiff(1:length(u),setdiff(Index,[1 length(u)]));
        Time4=Time3(Index3);
          
       
    for i1=CompoArtefact
         u(:,i1) = interp1(Time4,u(Index3,i1),Time3, 'linear');
    end
       
    Data(:,Index2)=transpose(u*s*v')+M;
            
    %Save as a newfile
    D=clone(D,['art_corr_' int_corr D.fname], [D.nchannels D.nsamples D.ntrials]);
    D(:,:,:)=Data;
    save(D);
end
