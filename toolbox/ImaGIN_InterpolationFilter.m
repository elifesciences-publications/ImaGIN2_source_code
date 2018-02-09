function D=ImaGIN_InterpolationFilter(S)

% Correct artefacts in frequency space
% Ref: Ramirez & Baillet, 2011

% S.method = 'linear' or 'spline' or 'other'

% try
    t=S.Fname;
% catch
%     t = spm_select(Inf, '\.mat$', 'Select data file');
% end
for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Time=time(D);
    
    
    F  = spm_figure('GetWin','Interactive');
    figure(F);clf
    
    
%     try
        EventType = S.EventType;
%     catch
%         EventType=spm_input('Events of interest', '+1', 'i',1);
%     end
%     try
        StartInterpolation = S.StartInterpolation;
%     catch
%         StartInterpolation=spm_input('Start of interpolation [sec]', '+1', 'r',-0.001);
%     end
%     try
        EndInterpolation = S.EndInterpolation;
%     catch
%         EndInterpolation=spm_input('End of interpolation [sec]', '+1', 'r',0.003);
%     end
    StartInterpolation=round(StartInterpolation*fsample(D));
    EndInterpolation=round(EndInterpolation*fsample(D));
    
    Bad=badchannels(D);
    Nchan=nchannels(D);
    Good=setdiff(1:Nchan,Bad);
    Nchan=length(Good);
    
    ev=events(D);
    
    Data=D(:,:);
    
    switch S.method
        
        case 'linear'
            
            for i1=1:length(ev)
                if mod(i1,200)==0
                    disp([i1 length(ev)])
                end
                %         if ~isempty(intersect(EventType,ev(i1).value))
                if ~isempty(intersect(EventType,ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    tmp1=D(:,ind+[StartInterpolation EndInterpolation]);
                    tmp2=zeros(size(Data,1),length([StartInterpolation:EndInterpolation]));
                    for i2=1:size(Data,1)
                        tmp2(i2,:)=interp1([StartInterpolation EndInterpolation],tmp1(i2,:),[StartInterpolation:EndInterpolation]);
                    end
                    Data(:,ind+[StartInterpolation:EndInterpolation])=tmp2;
                end
            end
            
            
        case 'spline'
            
            ind1=[];
            ind2=[];
            for i1=1:length(ev)
                %         if ~isempty(intersect(EventType,ev(i1).value))
                if ~isempty(intersect(EventType,ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    ind1=[ind1 ind+StartInterpolation];
                    ind2=[ind2 ind+EndInterpolation];
                end
            end
            ind1=sort(ind1);
            ind2=sort(ind2);
            Time=time(D);
            Index=1:length(Time);
            for i1=1:length(ind1)
                Index=setdiff(Index,ind1(i1):ind2(i1));
            end
            Time2=Time(Index);
            
            for i1=1:size(Data,1)
                Data(i1,:) = spline(Time2,Data(i1,Index),Time);
            end
            
            
            
        case 'other'
            
            NCompo=2;
            CCthresh=0.5;
            ind12=[];
            ind22=[];
            for i1=1:length(ev)
                %         if ~isempty(intersect(EventType,ev(i1).value))
                if strcmp(EventType,strvcat(ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    ind12=[ind12 ind+StartInterpolation*3];
                    ind22=[ind22 ind+EndInterpolation*3];
                end
            end
            ind12=sort(ind12);
            ind22=sort(ind22);
            Index2=[];
            for i1=1:length(ind12)
                Index2=[Index2 ind12(i1):ind22(i1)];
            end
            [u s v]=svd(Data(:,Index2)',0);
            
            f=abs(fft(u));
            CC=zeros(1,size(f,2));
            for i1=1:size(f,2)
                cc=corrcoef(f(:,1),f(:,i1));
                CC(i1)=cc(2);
            end
            CompoArtefact=find(CC>CCthresh);
            %     s=diag(s);
            %     s(CompoArtefact)=0;
            %     s=diag(s);
            
            
            U=u(:,CompoArtefact);
            windowSize = round(1*size(U,1)/length(ev));
            for i1=1:length(CompoArtefact)
                U2(:,i1)=2*U(:,i1)-filter(ones(1,windowSize)/windowSize,1,U(:,i1))-flipud(filter(ones(1,windowSize)/windowSize,1,flipud(U(:,i1))));
            end
            Artefact=ImaGIN_normalisation(sum(abs(U2),2),1);
            windowSize = round(1*length(Artefact)/length(ev));
            Artefact=2*Artefact-filter(ones(1,windowSize)/windowSize,1,Artefact)-flipud(filter(ones(1,windowSize)/windowSize,1,flipud(Artefact)));
            
            Index3=setdiff(1:length(Artefact),setdiff(find(Artefact>0.5),[1 length(Artefact)]));
            Time3=Time(Index2);
            Time4=Time3(Index3);
            %     for i1=1:NCompo
            for i1=CompoArtefact
                %         u(:,i1) = spline(Time4,u(Index3,i1),Time3);
                u(:,i1) = interp1(Time4,u(Index3,i1),Time3,'linear');
            end
            Data(:,Index2)=transpose(u*s*v');
            
            
    end
    
    
    %Save as a newfile
    D=clone(D,['i' D.fname], [D.nchannels D.nsamples D.ntrials]);
    D(:,:,:)=Data;
    save(D);
end
