function [Stimulation,StimulationIndex,StimulationFreqU]=ImaGIN_StimDetect(S)

try
    Filename=S.Fname;
catch
    Filename = spm_select(1, '\.mat$', 'Select data file');
end

D=spm_eeg_load(Filename);

try
    SelChan=S.Channels;
    if isempty(SelChan)
        SelChan=setdiff(setdiff(1:nchannels(D),indchantype(D,'ECG')),badchannels(D));
    end
catch
    SelChan = spm_input('Select channel', 1, 'i');
end

try
    Start=S.StimStart;
catch
    Start=[];
end

try
    End=S.StimEnd;
catch
    End=[];
end

try
    Stim=S.StimFreq;
catch
    Stim=spm_input('Stimulation frequency [Hz]', '+1', 'r');
end

try
    FindBadChannels=S.FindBadChannels;
catch
    FindBadChannels=1;
end

try
    StimContinuous=S.StimContinuous;
catch
    StimContinuous=spm_input('Continuous stimulation ', '+1', 'Yes|No');
    switch StimContinuous
        case 'Yes'
            StimContinuous=1;
        case 'No'
            StimContinuous=0;
    end
end

try
    FileOut=S.FileOut;
catch
    FileOut = [Filename(1:end-3) 'txt'];
end


% D=timeonset(D,0);
% save(D)
Time=time(D);


StimInit=Stim;
Stim=D.fsample/StimInit;
if isempty(Start)
    Start=1;
else
    Start=min(indsample(D,Start));
end
if isempty(End)
    End=nsamples(D);
else
    End=min(indsample(D,End));
    if isnan(End)
        End=nsamples(D);
    end
end

% Data=sum(D(SelChan,Start:End),1);
Data=D(SelChan,Start:End);

%find bad channels which saturate
if FindBadChannels
    L=zeros(1,size(Data,1));
    for i1=1:size(Data,1)
        tmp=abs(Data(i1,1:2:end));
        L(i1)=length(find(tmp==max(tmp)))/length(tmp);
    end
%     GoodChannels=find(L<200/size(Data,2));
    GoodChannels=find(L<0.05);      %Assume that bad channels saturate 5% of time
    Data=ImaGIN_normalisation(Data(GoodChannels,:),2);
end

%Remove line noise at 50 Hz
for i1=1:size(Data,1)
    Data(i1,:)=ImaGIN_notch(Data(i1,:),D.fsample,50,250);
end


%Simple version based on averaging between electrodes
stimulation = [];
if 1==1
    
    % d=[0 abs(diff(Data,2)) 0];
    d=[zeros(size(Data,1),1) abs(diff(Data,2,2)) zeros(size(Data,1),1)];
%     d=ImaGIN_normalisation(d,2);
%     [tmp,Order]=sort(max(d'),'descend');
%     d=mean(d(Order(1:ceil(length(Order)/2)),:));


    d=ImaGIN_normalisation(d,2);
    d1=max(d,[],1);
    d1=ImaGIN_normalisation(d1,2);
    d2=mean(d,1);
    d2=ImaGIN_normalisation(d2,2);
    
    d=d1+d2;

    
    d=ImaGIN_normalisation(d,2);
    
    % [tmp1,tmp2]=max(d);
    % Index=find(d>tmp1/4);
    Index=find(d>2);
    Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
    if ~isempty(Index)
        IndexClust=zeros(size(Index));
        IndexClust(1)=1;
        for i1=2:length(Index)
            tmp1=find(abs(Index(i1)-Index)<ceil(Stim/2));
            if max(IndexClust(tmp1))==0
                IndexClust(tmp1)=max(IndexClust(1:i1-1))+1;
            else
                IndexClust(tmp1)=max(IndexClust(tmp1));
            end
            %         IndexClust(tmp1(find(IndexClust(tmp1)==0)))=min(IndexClust(1:i1-1))+1;
        end
        IndexNew=zeros(1,max(IndexClust));
        for i1=1:length(IndexNew)
            tmp=find(IndexClust==i1);
            IndexNew(i1)=Index(tmp(floor(median(1:length(tmp)))));
        end
        Index=IndexNew;
        
        Template=0;
        for i1=1:length(Index)
            Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
        end
        
        % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
        
        
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        %Threshold cc
        Threshold=0.4;
        Threshold=min([0.4 median(cc(find(cc~=0)))]);
        tmp=find(cc>Threshold);
        stimulation=tmp;
        
        
        %Second pass
        Index=stimulation;
        Template=0;
        for i1=1:length(Index)
            Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
        end
        Index=find(d>2);
        Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        %Threshold cc
        Threshold=0.4;
        tmp=find(cc>Threshold);
        stimulation=tmp;

    end 
    
    if ~isempty(stimulation)
        %try to correct onset inpired from Lena's method
        tmps=sign(mean(Data(:,stimulation),2));
        Data2=Data;
        for i1=1:size(Data2,1)
            Data2(i1,:)=tmps(i1)*Data2(i1,:);
        end
        if 1==0
            %Same offset for all pulses
            TemplateData=0;
            for i1=1:length(stimulation)
                TemplateData=TemplateData+mean(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]),1);
            end
            TimeTemplate=[-ceil(Stim/2):ceil(Stim/2)]./fsample(D);
            TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.02&TimeTemplate>-0.1));
            tmpoffset=find(TemplateNorm>10&TimeTemplate>-0.02);
            if isempty(tmpoffset)
                offset=0;
            else
                i1=min(tmpoffset);
                offset=ceil(Stim/2)+1-i1;
            end
        else
            %different offset for all pulses
            TimeTemplate=[-ceil(Stim/2):ceil(Stim/2)]./fsample(D);
            offset=stimulation;
            for i1=1:length(stimulation)
                TemplateData=mean(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]),1);
                TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.02&TimeTemplate>-0.1));
                tmpoffset=find(TemplateNorm>5&TimeTemplate>-0.02);
                if isempty(tmpoffset)
                    offset(i1)=0;
                else
                    i2=min(tmpoffset);
                    offset(i1)=ceil(Stim/2)+1-i2;
                end
                if abs(offset(i1))>5
                    offset(i1)=0;
                end                    
            end
        end
        
        stimulation=stimulation-offset;
    end
    
else
    
    %more complex version where start and end stim are automatically found
    
    clear stimulation StartStim EndStim template
    n=0;
    for i0=1:size(Data,1)
        % d=[0 abs(diff(Data,2)) 0];
        d=[0 abs(diff(Data(i0,:),2,2)) 0];
        d=ImaGIN_normalisation(d,2);
        
        % [tmp1,tmp2]=max(d);
        % Index=find(d>tmp1/4);
        Index=find(d>3);
        
        Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
        Template=0;
        for i1=1:length(Index)
            Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
        end
        
        % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
        
        
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        
        %Threshold cc
        Threshold=0.25;
        tmp=find(cc>Threshold);
        stimulation{i0}=tmp;
        if ~isempty(tmp)
            StartStim(i0)=min(tmp);
            EndStim(i0)=max(tmp);
            n=n+1;
            template(n,:)=Template;
        end
        
    end
    
    Template=mean(template);
    CC=zeros(1,size(template,1));
    for i0=1:length(CC)
        cc=corrcoef(template(i0,:),Template);
        CC(i0)=cc(2);
    end
    Index=find(StartStim>0);
    Index=Index(find(CC>0.6));
    
        
% % %     StartStim=median(StartStim(find(EndStim>0)))-fsample(D);
% % %     EndStim=median(EndStim(find(EndStim>0)))+fsample(D);
% %     StartStim=min(StartStim(find(EndStim>0)))-fsample(D);
% %     EndStim=max(EndStim(find(EndStim>0)))+fsample(D);
%     
%     
%     Duration=(EndStim-StartStim)./fsample(D);
%     Index=find(StartStim>0);
% %     idx=kmeans(Duration(Index),2);
%     idx=kmeans([StartStim(Index);EndStim(Index);Duration(Index)]',2);
%     if mean(Duration(Index(find(idx==1))))>mean(Duration(Index(find(idx==2))))
%         Index=Index(find(idx==1));
%     else
%         Index=Index(find(idx==2));
%     end
%     
    
    
    StartStim=min(StartStim(Index))-fsample(D);
    EndStim=max(EndStim(Index))+fsample(D);
    
    
    d=[zeros(size(Data(Index,:),1),1) abs(diff(Data(Index,:),2,2)) zeros(size(Data(Index,:),1),1)];
    % d=ImaGIN_normalisation(d,2);
    % [tmp,Order]=sort(max(d'),'descend');
    % d=mean(d(Order(1:ceil(length(Order)/2)),:));
    d=mean(d);
    d=ImaGIN_normalisation(d,2);
    
    % [tmp1,tmp2]=max(d);
    % Index=find(d>tmp1/4);
    Index=find(d>3);
    
    Index=Index(find(Index>=StartStim&Index<=EndStim));
    
    Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
    Template=0;
    for i1=1:length(Index)
        Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
    end
    
    % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
    
    
    cc=zeros(size(d));
    % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
    for i1=Index
        tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
        cc(i1)=tmp(2);
        tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
        if length(tmp1)>1
            cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
        end
    end
    
    %Threshold cc
    Threshold=0.5;
    tmp=find(cc>Threshold);
    stimulation=tmp;
        
end




%remove outliers close to other stims
if length(stimulation)>1
    remove=[];
    for i1=1:length(stimulation)
        difftmp=abs(stimulation-stimulation(i1));
        if isempty(find(difftmp>0.95*Stim & difftmp<1.05*Stim))
            remove=[remove i1];
        end
    end
%     StimulationFreqU=D.fsample/median(diff(stimulation));   %uncorrected stim frequency
    stimulation=stimulation(setdiff(1:length(stimulation),remove));
end

if isempty(stimulation)
    StimulationFreqU=[];
    StimulationIndex = [];
    Stimulation = [];
else
    Index=find(d(stimulation(1):stimulation(end))>2);
    Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
    if length(Index>1) > 2 % Index should be great than 2 not 0
        tmp=diff(Index);
        tmp=tmp(find(tmp>max([2 D.fsample/100])));
        if ~isempty(tmp)
            tmp=sort(tmp);
            tmp=tmp(ceil(length(tmp)/3):end);
            StimulationFreqU=D.fsample/median(tmp);   %uncorrected stim frequency
        else
            StimulationFreqU = 0;
        end
    else
        StimulationFreqU=0;   %uncorrected stim frequency
    end
    
    
    
    %fill the gaps
    if StimContinuous
        tmp=stimulation;
        stimulation=tmp(1);
        while stimulation(end)<tmp(end)
            [tmp1 tmp2]=min(abs(tmp-(stimulation(end)+Stim)));
            if tmp1>Stim/2
                stimulation=[stimulation stimulation(end)+Stim];
            else
                stimulation=[stimulation tmp(tmp2)];
            end
        end
    end
    
    StimulationIndex=stimulation+Start-1;
    Stimulation=Time(round(StimulationIndex));
end

% 
% tmp2=diff(tmp);
% tmp3=[1 find(tmp2~=1)+1 length(tmp)];
% Stimulation=zeros(length(tmp3)-1,1);
% for i1=1:length(tmp3)-1
%     [t1,t2]=max(cc(tmp(tmp3(i1):tmp3(i1+1))));
%     Stimulation(i1)=tmp(t2+tmp3(i1)-1);
% end
% Stimulation=Time(Start-1+Stimulation);
% 
% ok=1;
% while ok
%     StimulationOK=ones(size(Stimulation));
%     for i1=2:length(Stimulation)
%         if Stimulation(i1)-Stimulation(i1-1)<1/(10*StimInit)
%             StimulationOK(i1)=0;
%         end
%     end
%     Index=find(StimulationOK==0);
%     if isempty(Index)
%         ok=0;
%     else
%         Stimulation=Stimulation(find(StimulationOK));
%     end
% end


% P = spm_str_manip(Filename, 'H');
% F = spm_str_manip(Filename, 'tr');

fid=fopen(FileOut,'w');
fprintf(fid,'%f\n',Stimulation);
fclose(fid);

% Filename=fullfile(DirOut,[F '_StimulationIndex.txt']);
% fid=fopen(Filename,'w');
% fprintf(fid,'%f\n',StimulationIndex);
% fclose(fid);


