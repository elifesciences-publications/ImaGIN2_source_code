function ImaGIN_StimCorrec(S)

try
    Filename=S.Filename;
catch
    Filename = spm_select(1, '\.mat$', 'Select data file');
end
D=spm_eeg_load(Filename);
if isfield(D,'time')
    time=D.time;
else
    time=[0:1/D.fsample:(D.nsamples-1)/D.fsample];
    time=time-D.timeonset;
end

try
    SelChan=S.Channels;
catch
    SelChan = spm_input('Select channel', 1, 'i',sprintf('1:%d',D.nchannels));
end

try
    Start=S.StimStart;
catch
    Start=spm_input('Start of analysis window [sec]', '+1', 'r',min(time));
end

try
    End=S.StimEnd;
catch
    End=spm_input('End of analysis window [sec]', '+1', 'r',max(time));
end

try
    Stim=S.StimFreq;
catch
    Stim=spm_input('Stimulation frequency [Hz]', '+1', 'r');
end

try
    FlagContinuous=S.Continuous;
catch
    FlagContinuous=spm_input('Continuous stimulation','+1','Yes|No',[1 0]);
end


if isempty(Start)
    Start=1;
else
    Start=unique(find(abs(time-Start)==min(abs(time-Start))));
end
if isempty(End)
    End=D.nsamples;
else
    End=unique(find(abs(time-End)==min(abs(time-End))));
end
StimInit=Stim;
Wo = StimInit*(time(2)-time(1))*2;  BW = Wo/35;
[bnotch,anotch] = iirnotch(Wo,BW);  

Data=D(SelChan,Start:End);
data=abs(sum(Data)-filter(ones(1,2)/2,1,sum(Data)));
Time=D.time(Start:End);
if FlagContinuous
    Threshold=1*mean(data);
else
    Threshold=3*mean(data);
end
[SpikeRate,MinLoc,MaxLoc]=ImaGIN_SpikeDetect(data,Time,StimInit,Threshold,[],0,FlagContinuous);

%use svd
ISI=median(diff(MaxLoc));
Gap=min([round(ISI/4) round(0.05*D.fsample)]);
Index=MaxLoc'*ones(1,ISI)+ones(length(MaxLoc),1)*[-Gap:ISI-Gap-1];
DataCorrected=Data;
for i1=1:size(Data,1)
    tmp=Data(i1,:);
    tmp=tmp(Index);
    M=mean(tmp,2)*ones(1,size(tmp,2));
    tmp=tmp-M;
    [u,s,v]=svd(tmp,0);
    s=diag(s);
    %Zeroing components
% %     variance=cumsum(s.^2);variance=variance/max(variance);
% %     Ncompo=length(find(variance<0.98));
% %     Ncompo=round(ISI/2);
% %     s(1:Ncompo)=0;
%     for i2=1:size(v,2)
%         c=v(:,i2);
%         [tmp1,tmp2]=max(abs(c));
%         if tmp2<=2*Gap+1&tmp2>1
%             s(i2)=0;
%         end
%     end
    %Removing artefact by linear interpolation
    for i2=1:size(v,2)
        c=v(:,i2);
        [tmp1,tmp2]=max(abs(detrend(c)));
        if tmp2<=2*Gap+1&tmp2>1
            v(:,i2) = interp1([1 2*Gap+2:length(c)],c([1 2*Gap+2:length(c)]),1:length(c),'linear');
        end
    end
    correct=u*diag(s)*v';
    correct=correct+M;
    DataCorrected(i1,Index(:))=correct(:);
    %smooth
    DataCorrected(i1,:)=filter(ones(1,Gap)/Gap,1,DataCorrected(i1,:));
    DataCorrected(i1,:)=filter(bnotch,anotch,DataCorrected(i1,:));
end

%Save as a newfile
Data=D(:,:);
Data(SelChan,Start:End)=DataCorrected;
D.data=[];
P = spm_str_manip(Filename, 'H');
Prefix='c';
D.fname=[Prefix '_' D.fname];
D.fnamedat=[Prefix '_' D.fnamedat];
fpd = fopen(fullfile(P, D.fnamedat), 'w');
D.scale=zeros(D.Nchannels,1,D.Nevents);
for i1=1:D.Nevents
    D.scale(:,1,i1) = spm_eeg_write(fpd, Data(:,:,i1), 2, D.datatype);
end
fclose(fpd);
if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end
