function ImaGIN_SeizureDetect(S)
% Compute instantaneous power with time-windowed fft to detect seizures
% 
% USAGE:  D = ImaGIN_SeizureDetect(S)

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

try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select EEG mat file');
end
D=spm_eeg_load(DD);
P=spm_str_manip(DD,'h');


try
    Job=S.Job;
catch
    Ctype = {
        'Create',...
        'Events',...
        'Statistics',...
        'Delete'};
    str   = 'Seizure manipulation ';
    Sel   = spm_input(str, 1, 'm', Ctype);
    Job = Ctype{Sel};
end
switch Job
    case 'Create'
        try
            ImaGIN_SeizureSeizureCreate(D,P,S);
        catch
            ImaGIN_SeizureSeizureCreate(D,P);
        end
    case 'Replace'
        try
            ImaGIN_SeizureSeizureReplace(D,P,S);
        catch
            ImaGIN_SeizureSeizureReplace(D,P);
        end
    case 'Events'
       try
            ImaGIN_SeizureSeizureEvents(D,P,S);
       catch
           ImaGIN_SeizureSeizureEvents(D,P);
       end
    case 'Statistics'
        try
            ImaGIN_SeizureSeizureStatistics(D,P,S);
        catch
            ImaGIN_SeizureSeizureStatistics(D,P);
        end
    case 'Delete'
        try
            ImaGIN_SeizureSeizureDelete(D,P,S);
        catch
            ImaGIN_SeizureSeizureDelete(D,P);
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_SeizureSeizureCreate(D,P,S)

if isfield(D,'Seizure')   %Recalculate seizure
    if ~isempty(D.Seizure)
        KeepNumber=length(D.Seizure)+1;
        NameNumber=str2num(D.Seizure{end}.Name(end))+1;
    else
        KeepNumber=1;
        NameNumber=1;
    end
else
    D.Seizure=[];
    KeepNumber=1;
    NameNumber=1;
end

Time=time(D);

try
    SelChan = S.Channel;
catch
    if isfield(D,'Seizure')
        SelChan = spm_input('Select channel(s) ', '+1', 'i');
    else
        SelChan = spm_input('Select channel(s) ', 1, 'i');
    end
end

try
    Method = S.Method;
catch
    Method = spm_input('Measure for seizure detection ',1,'SVD|Spike|DC|Entropy|Amplitude');
end

switch Method
    case 'SVD'
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r',[num2str(min(Time)) ':' num2str(max(Time))]);
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2);
        end
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
       
        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
    
    case 'Entropy'
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r',[num2str(min(Time)) ':' num2str(max(Time))]);
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2);
        end
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
        try
            EmbeddingDimension=S.EmbeddingDimension;
        catch
            EmbeddingDimension = spm_input('Embedding dimension', '+1', 'i',3);
        end
        try
            Subject=S.Subject;
        catch
            Subject = spm_input('Subject', '+1', 'GAERS|None');
        end
       
        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
        D.Seizure{KeepNumber}.EmbeddingDimension=EmbeddingDimension;

    case 'DC'
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r',[num2str(min(Time)) ':' num2str(max(Time))]);
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2);
        end
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;

        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;

    case 'Spike'
        try
            Start = S.Start;
        catch
            Start=spm_input('Start of analysis window [sec]', '+1', 'r',Time(1));
        end
        try
            End = S.End;
        catch
            End=spm_input('End of analysis window [sec]', '+1', 'r',Time(end));
        end
        
        try
            Freq=S.frequencies;
        catch
            Freq = spm_input('Frequency of interest [Hz]', '+1', 'r', '', 1);
        end

        try
            ThreshCC=S.ThreshCC;
        catch
            ThreshCC = spm_input('Threshold on spike correlation (>0, 0=auto) ', '+1', 'r',0,1);
        end

        try
            ThreshData=S.ThreshData;
        catch
            ThreshData = spm_input('Threshold on amplitude (>0, 0=auto) ', '+1', 'r',0,1);
        end

        try
            Coarse=S.Coarse;
        catch
            Coarse = spm_input('Downsampling ', '+1', 'i', 0);
        end
        
        if ThreshData<=0
            try
                Baseline=S.Baseline;
            catch
                Baseline = spm_input('Baseline time window [s]', '+1', 'r', '', 2);
            end
            if isempty(Baseline)
                ThreshData=std(D(SelChan,:));
            else
                tmp=sort(abs(D(SelChan,find(Time>=Baseline(1),1):find(Time<=Baseline(2),1,'last'))));
                ThreshData=tmp(round(0.99*length(tmp)));
            end
            %Save parameter analysis
            D.Seizure{KeepNumber}.Baseline=Baseline;
        end
        %Save parameter analysis
        D.Seizure{KeepNumber}.ThreshData=ThreshData;
        D.Seizure{KeepNumber}.Freq=Freq;

    case 'Amplitude'
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r',[num2str(min(Time)) ':' num2str(max(Time))]);
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2);
        end
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
        try
            SpikeWidth=S.SpikeWidth;
        catch
            SpikeWidth = spm_input('Spike width [sec]', '+1', 'r',0.2);
        end

        try
            Baseline=S.Baseline;
        catch
            Baseline = str2num(spm_input('Baseline time window [s]', '+1', 's', ''));
        end
        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
        D.Seizure{KeepNumber}.SpikeWidth=SpikeWidth;
end


%Save parameter analysis
D.Seizure{KeepNumber}.Start=Start;
D.Seizure{KeepNumber}.End=End;
D.Seizure{KeepNumber}.SelChan=SelChan;
D.Seizure{KeepNumber}.Name=['Seizure ' num2str(NameNumber)];

%Crop the data
Start=unique(find(abs(Time-Start)==min(abs(Time-Start))));
End=unique(find(abs(Time-End)==min(abs(Time-End))));
Data=D(SelChan,Start:End);
Data(find(isnan(Data)))=0;
Time=Time(Start:End);


if exist('TimeWindow','var')
    TimeWindow=TimeWindow(find(TimeWindow-TimeWindowWidth/2>=Time(1) & TimeWindow+TimeWindowWidth/2<=Time(end)));
end

switch Method
    case 'SVD'
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            
            data=Data(:,win).*(ones(size(Data,1),1)*hanning(length(win))');
            data=ImaGIN_normalisation(data,2,[]);

            for i2=1:size(data,1)
                data(i2,:)=ImaGIN_bandpass(data(i2,:),fsample(D),7,90);
            end
            data=ImaGIN_normalisation(data,2,[]);
            [u,s,v]=svd(data',0);
            s=diag(s);
            
            Seizure(i1)=sum(s(1:ceil(length(s)/2)).^2)/sum(s.^2);
        end

        S = interp1(TimeWindow,Seizure,Time,'linear');
    
    
    case 'DC'        %search for DC shift (asymetry) during seizure
        %notch filter
        Wo = 50*(Time(2)-Time(1))*2;  
        BW = Wo/35;
        Data = ImaGIN_notch(Data, Wo, BW);

        Seizure=zeros(1,length(TimeWindow));
        if median(Data)<0
            PeakPos=1;  %detect positive peaks
        else
            PeakPos=0;
        end
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            
            data=Data(win);
            data=data-median(data);

            dmin2=sort(data);
            dmax2=sort(data,'descend');

            data=sort(data);
            data=data(round(0.3*length(data)):round(0.7*length(data)));

            if PeakPos
                Seizure(i1)=(sum(dmax2(find(dmax2>-dmin2(ceil(0.1*length(dmin2))))))+sum(dmin2(find(dmin2<dmin2(ceil(0.1*length(dmin2)))))))/std(data);
            else
                Seizure(i1)=-(sum(dmax2(find(dmax2>dmax2(ceil(0.1*length(dmax2))))))+sum(dmin2(find(dmin2<-dmax2(ceil(0.1*length(dmax2)))))))/std(data);
            end
        end

        tmp=sort(abs(Seizure(find(Seizure<0))));
        tmp=tmp(round(0.95*length(tmp)));
        Seizure=abs(Seizure)./tmp;

        S = interp1(TimeWindow,Seizure,Time,'linear');

    case 'Amplitude'        %search for spike amplitude
        %notch filter
        Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
        Data = ImaGIN_notch(Data, Wo, BW);

        if isempty(Baseline)
            Data=ImaGIN_normalisation(Data,2,[]);
        else
            Bsl=[];
            for i=1:length(Baseline)
                Bsl(i)=indsample(D,Baseline(i));
            end
            Data=ImaGIN_normalisation(Data,2,Bsl(1):Bsl(2));
        end
        SpikeWidthSamples=SpikeWidth*D.fsample;
        
        %find local Maxima & Minima
        MinLoc=zeros(size(Data));
        MaxLoc=zeros(size(Data));
        for i1=3:length(Data)-2
            if (Data(i1)==Data(i1-1)) && (Data(i1)==Data(i1+1))
                if (Data(i1)-Data(i1-2)<=0) && (Data(i1)-Data(i1+2)<=0)
                    MinLoc(i1)=1;
                end
                if (Data(i1)-Data(i1-2)>=0) && (Data(i1)-Data(i1+2)>=0)
                    MaxLoc(i1)=1;
                end
            else
                if (Data(i1)-Data(i1-1)<=0) && (Data(i1)-Data(i1+1)<=0)
                    MinLoc(i1)=1;
                end
                if (Data(i1)-Data(i1-1)>=0) && (Data(i1)-Data(i1+1)>=0)
                    MaxLoc(i1)=1;
                end
            end
        end
        MinLoc=find(MinLoc);
        MaxLoc=find(MaxLoc);
        
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            
            IndexMax=MaxLoc(find(MaxLoc>=win(1)&MaxLoc<=win(end)));
            IndexMin=MinLoc(find(MinLoc>=win(1)&MinLoc<=win(end)));
            %effets de bord
            if IndexMin(2)<IndexMax(1)
                IndexMin=IndexMin(2:end);
            end
            if IndexMax(2)<IndexMin(1)
                IndexMax=IndexMax(2:end);
            end
            
            for i2=1:length(IndexMax)
                if Data(IndexMax(i2))>0
                    tmp=find(abs(IndexMax(i2)-IndexMin)<SpikeWidthSamples&Data(IndexMin)<0);
                    S2=max(Data(IndexMax(i2))-Data(IndexMin(tmp)));
                    if S2>Seizure(i1)
                        Seizure(i1)=S2;
                    end
                end
            end
        end

        S = interp1(TimeWindow,Seizure,Time,'linear');
    
    case 'Entropy'        %permutation entropy
        %notch filter
        Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
        Data = ImaGIN_notch(Data, Wo, BW);

        Coarse=S.Coarse;
        TimeDelay=1;
        Seizure=ones(length(TimeDelay),length(TimeWindow));
        Power=ones(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            t=Time(win);
            f = ImaGIN_time2freq(t);
            switch  Subject
                case 'GAERS'
                    FreqInterest = find(f>=7&f<=11);
                    FreqNoInterest = find(f>=1&f<=6);
                case 'None'
                    FreqInterest = [];
                    FreqNoInterest = [];
            end
            power=abs(fft(Data(win)));
            Power(i1)=sum(Data(win).*Data(win));
            if isempty(FreqInterest)
                win=win(1:Coarse:floor(TimeWindowWidth*D.fsample));
                Seizure(i1)=ImaGIN_permutation_entropy(Data(win),EmbeddingDimension,TimeDelay);
            elseif mean(power(FreqInterest))>mean(power(FreqNoInterest))
                win=win(1:Coarse:floor(TimeWindowWidth*D.fsample));
                Seizure(i1)=ImaGIN_permutation_entropy(Data(win),EmbeddingDimension,TimeDelay);
            end
            if isnan(Seizure(i1))
            end
        end

        S = interp1(TimeWindow,Seizure,Time,'linear');

    case 'Spike'
        %notch filter
        Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
        Data = ImaGIN_notch(Data, Wo, BW);

        [S,Event1,Event2]=ImaGIN_spike_detect(Data,Time,Freq,ThreshData,ThreshCC,0,Coarse);
        D.Seizure{KeepNumber}.Event1=Event1;
        D.Seizure{KeepNumber}.Event2=Event2;
end

if Start>1
    S=[zeros(1,Start-1) S];
end
if End<D.nsamples
    S=[S zeros(1,D.nsamples-End)];
end

D.Seizure{KeepNumber}.data=S;

save(D);
spm('Pointer', 'Arrow');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_SeizureSeizureReplace(D,P,S)

NSeizure=length(D.Seizure);
if NSeizure>1
    str='';
    for i1=1:NSeizure
        str=[str num2str(i1) '|'];
    end
    str=str(1:end-1);
    KeepNumber = str2num(spm_input('Replace seizure time series number ',1,str));
else
    KeepNumber=1;
end

if isfield(D,'time')
    time=D.time;
else
    time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
    time=time-timeonset(D);
end

try
    SelChan = S.Channel;
catch
    if isfield(D,'Seizure')
        SelChan = spm_input('Select channel(s) ', '+1', 'i');
    else
        SelChan = spm_input('Select channel(s) ', 1, 'i');
    end
end

try
    Start = S.Start;
catch
    Start=spm_input('Start of analysis window [sec]', '+1', 'r',time(1));
end
try
    End = S.End;
catch
    End=spm_input('End of analysis window [sec]', '+1', 'r',time(end));
end

try
    Method = S.Method;
catch
    Method = spm_input('Measure for seizure detection ',1,'Spike|Power');
end

switch Method
    case 'Spike'
        try
            Freq=S.frequencies;
        catch
            Freq = spm_input('Frequency of interest [Hz]', '+1', 'r', '', 1);
        end

        try
            ThreshData=S.ThreshData;
        catch
            ThreshData = spm_input('Threshold on amplitude (>0, 0=auto) ', '+1', 'r',0,1);
        end

        if ThreshData<=0
            try
                Baseline=S.Baseline;
            catch
                Baseline = spm_input('Baseline time window [s]', '+1', 'r', '', 2);
            end
            if isempty(Baseline)
                ThreshData=std(D(SelChan,:));
            else
                tmp=sort(abs(D(SelChan,find(time>=Baseline(1),1):find(time<=Baseline(2),1,'last'))));
                ThreshData=tmp(round(0.99*length(tmp)));
            end
            %Save parameter analysis
            D.Seizure{KeepNumber}.Baseline=Baseline;
        end
        %Save parameter analysis
        D.Seizure{KeepNumber}.ThreshData=ThreshData;
end


%Save parameter analysis
D.Seizure{KeepNumber}.Start=Start;
D.Seizure{KeepNumber}.End=End;
D.Seizure{KeepNumber}.SelChan=SelChan;
D.Seizure{KeepNumber}.Freq=Freq;

try
    Bin=S.Bin;
catch
    Flag = spm_input('Detect start and end of seizures ','+1','Yes|No');
    switch Flag
        case 'Yes'
            Bin=1;
        case 'No'
            Bin=0;
    end
end

if Bin
    try
        SeizureInterval=S.SeizureInterval;
    catch
        SeizureInterval = spm_input('Minimum seizure interval [sec]', '+1', 'r',2,1);
    end
    try
        SeizureDuration=S.SeizureDuration;
    catch
        SeizureDuration = spm_input('Minimum seizure duration [sec]', '+1', 'r',2,1);
    end
    try
        ThreshDetect=S.ThreshDetect;
    catch
        ThreshDetect = spm_input('Threshold detection (>0, 0=auto) ', '+1', 'r',3.5,1);
    end
end

%Save parameter analysis
D.Seizure{KeepNumber}.Bin=Bin;
if Bin
    D.Seizure{KeepNumber}.SeizureInterval=SeizureInterval;
    D.Seizure{KeepNumber}.SeizureDuration=SeizureDuration;
    D.Seizure{KeepNumber}.ThreshDetect=ThreshDetect;
end

%Crop the data
Start=unique(find(abs(time-Start)==min(abs(time-Start))));
End=unique(find(abs(time-End)==min(abs(time-End))));
Data=D(SelChan,Start:End);
Time=time(Start:End);

switch Method
    case 'Spike'
        S=ImaGIN_spike_detect(Data,Time,Freq,ThreshData,0.5,1);
end

if Start>1
    S=[zeros(1,Start-1) S];
end
if End<D.nsamples
    S=[S zeros(1,D.nsamples-End)];
end

D.Seizure{KeepNumber}.data=S;

if Bin
    if ThreshDetect==0
        %Determine the threshold with a 4 Gaussian mixture model
        Coarse=max([1 floor(length(S)/2e3)]);
        X=ImaGIN_normalisation(S(1:Coarse:end),2)';
        [W,M,R,Tlogl] = gmmbvl_em(X,4,4,0,0,0);
        [M,order]=sort(M);
        W=W(order);
        R=R(order);
        x = linspace(min(X) - 3*max(R), max(X) + 3*max(R), 500 )';
        L = gmmbvl_em_gauss(x,M,R);
        L=L.*repmat(W',size(L,1),1);
        Thresh=x(max(find(L(:,end-1)>L(:,end))))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
    else
        Thresh=ThreshDetect;
    end

    D.Seizure{KeepNumber}.Thresh=Thresh;

    try
        tmp=find(S>=Thresh);
    catch
        tmp1=find(S(1:round(length(S)/3))>=Thresh);
        tmp2=round(length(S)/3)+find(S(round(length(S)/3)+[1:round(length(S)/3)])>=Thresh);
        tmp3=2*round(length(S)/3)+find(S(2*round(length(S)/3)+1:end)>=Thresh);
        tmp=[tmp1 tmp2 tmp3];
    end
    End=[tmp(find(diff(tmp)>1)) tmp(end)];
    Start=tmp([1 find(diff(tmp)>1)+1]);
    SeizureInterval=SeizureInterval*D.fsample;
    SeizureDuration=SeizureDuration*D.fsample;

    %Remove very small bits (<0.2 s or 2 samples)
    Remove=[];
    for i1=1:length(Start)
        try
            if End(i1)-Start(i1)<max([0.2*D.fsample 2])
                Remove=[Remove i1];
            end
        end
    end
    Start=Start(setdiff(1:length(Start),Remove));
    End=End(setdiff(1:length(End),Remove));

    RemoveStart=[];
    RemoveEnd=[];
    for i1=1:length(Start)
        try
            if Start(i1+1)-End(i1)<SeizureInterval
                RemoveStart=[RemoveStart i1+1];
                RemoveEnd=[RemoveEnd i1];
            end
        end
    end
    Start=Start(setdiff(1:length(Start),RemoveStart));
    End=End(setdiff(1:length(End),RemoveEnd));
    Remove=[];
    for i1=1:length(Start)
        try
            if End(i1)-Start(i1)<SeizureDuration
                Remove=[Remove i1];
            end
        end
    end
    Start=Start(setdiff(1:length(Start),Remove));
    End=End(setdiff(1:length(End),Remove));
    if length(Start)>length(End)
        Start=Start(1:end-1);
    end

    %add events
    Events=events(D);
    if(~isempty(Events))
        try
            Events(1).type;
        catch
            Events=Events{1};
        end
    end
    Index=D.Seizure{KeepNumber}.Name(end);
    NeventOld=size(Events,2);
    NewName{1}=['Seizure' Index 'Start'];
    NewName{2}=['Seizure' Index 'End'];
    Timing{1}=Start;
    Timing{2}=End;
    for i1=1:2
        Events(i1+NeventOld).type=NewName{i1};
        Events(i1+NeventOld).value=i1+NeventOld;
        Events(i1+NeventOld).time=time(Timing{i1});
    end
    Events(1+NeventOld).duration=time(Timing{i1})-time(Timing{i1});

    D=events(D,1,Events);
end
save(D);
spm('Pointer', 'Arrow');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_SeizureSeizureStatistics(D,P,S)

NSeizure=length(D.Seizure);
if NSeizure>1
    str='';
    for i1=1:NSeizure
        str=[str num2str(i1) '|'];
    end
    str=str(1:end-1);
    KeepNumber = str2num(spm_input('Write statistics of seizure number ',1,str));
else
    KeepNumber = 1;
end

Index=D.Seizure{KeepNumber}.Name(end);
NewName{1}=['Seizure' Index 'Start'];
NewName{2}=['Seizure' Index 'End'];
Name{1}='_Seizure';
if D.Seizure{KeepNumber}.Somnolence
    NewName2{1}=['Somnolence' Index 'Start'];
    NewName2{2}=['Somnolence' Index 'End'];
    Name{2}='_Somnolence';
end

for i=1:size(Name,2)
    Data=[];
    Events=events(D);
    if i==2
        NewName=NewName2;
    end
    for i1=1:2
        Data1=[];
        for i2=1:size(Events,2)
            if strcmp(Events(i2).type,NewName{i1})
                Data1=[Data1;Events(i2).time];
            end
        end
        Data=[Data Data1];
    end
    Data=[Data Data(:,2)-Data(:,1)];

    %Write text file
    name=fname(D);
    File=fullfile(D.path,[name(1:end-4) Name{i} Index '.txt']);
    fid = fopen(File,'wt');
    fprintf(fid,' Start [s] /  End [s] / Duration [s] / Mean duration [s] / Std duration [s] / Total duration [s] / Rate [/h]\n');
    fprintf(fid,'\n');
    fprintf(fid,[' ' convertStoHHMMSS(Data(1,1)) '    ' convertStoHHMMSS(Data(1,2)) '    ' convertStoHHMMSS(Data(1,2)-Data(1,1))  '         '  convertStoHHMMSS(mean(Data(:,2)-Data(:,1))) '           ' convertStoHHMMSS(std(Data(:,2)-Data(:,1))) '            ' convertStoHHMMSS(sum(Data(:,2)-Data(:,1))) '            ' num2str(size(Data,1)/(D.nsamples/(D.fsample*3600))) '\n']);
    for i1=2:size(Data,1)
        fprintf(fid, [' ' convertStoHHMMSS(Data(i1,1)) '    ' convertStoHHMMSS(Data(i1,2)) '    ' convertStoHHMMSS(Data(i1,2)-Data(i1,1)) '\n']);
    end
    fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function convertStoHHMMSS(t1)
datestr(t1/(24*60*60), 'HH:MM:SS');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_SeizureSeizureEvents(D,P,S)

NSeizure=length(D.Seizure);
CombineFlag=0;
if NSeizure>1
    CombineFlag = spm_input('Combine two measures ',1,'Yes|No');
    if strcmp(CombineFlag,'Yes')
        CombineFlag=1;
    else
        CombineFlag=0;
    end
   
    
    str='';
    for i1=1:NSeizure
        str=[str num2str(i1) '|'];
    end
    str=str(1:end-1);
    KeepNumber = str2num(spm_input('Put events on seizure number ','+1',str));
    
    if CombineFlag
        KeepNumber2 = str2num(spm_input('Using seizure number ','+1',str));
    end
        
else
    KeepNumber = 1;
end

Time=time(D);

try
    SeizureInterval=S.SeizureInterval;
catch
    SeizureInterval = spm_input('Minimum seizure interval [sec]', '+1', 'r',1,1);
end
try
    SeizureDuration=S.SeizureDuration;
catch
    SeizureDuration = spm_input('Minimum seizure duration [sec]', '+1', 'r',5,1);
end
try
    ThreshDetect=S.ThreshDetect;
catch
    ThreshDetect = spm_input('Low threshold detection (>0, 0=auto) ', '+1', 'r',0,1);
end
try
    ThreshDetectAmp=S.ThreshDetectAmp;
catch
    ThreshDetectAmp = spm_input('High threshold detection (Amp, >0, 0=no) ', '+1', 'r',0,1);
end
try
    UseFrequency=S.UseFrequency;
catch
    UseFrequency=spm_input('Use frequencies ',1,'Yes|No');
    if strcmp(UseFrequency,'Yes')
        UseFrequency=1;
    else
        UseFrequency=0;
    end
end

try
    Somnolence=S.Somnolence;
catch
    Somnolence=spm_input('Detect somnolences ',1,'Yes|No');
    if strcmp(Somnolence,'Yes')
        Somnolence=1;
    else
        Somnolence=0;
    end
end
D.Seizure{KeepNumber}.Somnolence=Somnolence;

if CombineFlag
    try
        ThreshDetect2=S.ThreshDetect2;
    catch
        ThreshDetect2 = spm_input('Threshold for 2nd measure ', '+1', 'r',6,1);
    end
end

CombineAandDC=0;
if NSeizure>1
    CombineAandDC = spm_input('Combine amplitude and DC measures ',1,'Yes|No');
    if strcmp(CombineAandDC,'Yes')
        CombineAandDC=1;
        NumberDC = str2num(spm_input('DC on seizure number ','+1',str));
        ThreshDetect2 = spm_input('Threshold for DC measure ', '+1', 'r',1,1);
    else
        CombineAandDC=0;
    end
end

S=D.Seizure{KeepNumber}.data;

if ThreshDetect==0
    SS=S(find(~isnan(S)))';

    %use kmeans assuming 2 clusters
    y=ImaGIN_kMeansCluster(SS,2);
    c=zeros(1,2);
    c(1)=mean(y(find(y(:,end)==1)));
    c(2)=mean(y(find(y(:,end)==2)));
    Thresh = min(c)+std(c);
else
    Thresh=ThreshDetect;
end
ThreshAmp=ThreshDetectAmp;

%Save parameter analysis
D.Seizure{KeepNumber}.Bin=1;
D.Seizure{KeepNumber}.SeizureInterval=SeizureInterval;
D.Seizure{KeepNumber}.SeizureDuration=SeizureDuration;
D.Seizure{KeepNumber}.ThreshDetect=Thresh;
D.Seizure{KeepNumber}.ThreshDetectAmp=ThreshAmp;
if CombineFlag
    D.Seizure{KeepNumber}.ThreshDetect2=ThreshDetect2;
    S2=D.Seizure{KeepNumber2}.data;
end

tmp=find(S>=Thresh);
End=[tmp(find(diff(tmp)>1)) tmp(end)];
Start=tmp([1 find(diff(tmp)>1)+1]);
SeizureInterval=SeizureInterval*D.fsample;
SeizureDuration=SeizureDuration*D.fsample;

N=1;
while N~=0
    %Classify bits (assume large bits >MinSize s or 2 samples
    MinSize=2;
    Flag=zeros(1,length(Start));
    for i1=1:length(Start)
        try
            if End(i1)-Start(i1)<max([MinSize*D.fsample 2])
                Flag(i1)=1; %small bits
            end
        end
    end
    %Agregate small bits (<MinSize s or 2 samples) if the nearest large bits is closer than 1 s
    Th=1;
    RemoveStart=[];
    RemoveEnd=[];
    Flag0=find(Flag==0);
    for i1=1:length(Start)
        if Flag(i1)
            [tmp,tmp1]=min(abs(Start(i1)-End(Flag0))/D.fsample);
            if tmp>Th
                [tmp,tmp1]=min(abs(End(i1)-Start(Flag0))/D.fsample);
                if tmp<Th
                    RemoveStart=[RemoveStart i1+1:Flag0(tmp1)];
                    RemoveEnd=[RemoveEnd i1:Flag0(tmp1)-1];
                end
            else 
                RemoveStart=[RemoveStart Flag0(tmp1)+1:i1];
                RemoveEnd=[RemoveEnd Flag0(tmp1):i1-1];
            end
        end
    end
    N=length(RemoveStart);
    Start=Start(setdiff(1:length(Start),RemoveStart));
    End=End(setdiff(1:length(End),RemoveEnd));
end

%Classify bits (assume large bits >MinSize s or 2 samples
Flag=zeros(1,length(Start));
for i1=1:length(Start)
    try
        if End(i1)-Start(i1)<max([MinSize*D.fsample 2])
            Flag(i1)=1; %small bits
        end
    end
end
%Remove very small bits (<MinSize s or 2 samples) if the nearest large bits is further than 1 s
Remove=find(Flag);
Start=Start(setdiff(1:length(Start),Remove));
End=End(setdiff(1:length(End),Remove));

RemoveStart=[];
RemoveEnd=[];
for i1=1:length(Start)
    try
        if Start(i1+1)-End(i1)<SeizureInterval
            RemoveStart=[RemoveStart i1+1];
            RemoveEnd=[RemoveEnd i1];
        end
    end
end
Start=Start(setdiff(1:length(Start),RemoveStart));
End=End(setdiff(1:length(End),RemoveEnd));
Remove=[];
if ThreshAmp==0 %only test duration
    try
        if End(i1)-Start(i1)<SeizureDuration
            Remove=[Remove i1];
        end
    end
else            %test duration and amplitude
    for i1=1:length(Start)
        try
            if End(i1)-Start(i1)<SeizureDuration || ~isempty(find(S(Start(i1):End(i1))>=ThreshAmp))
                Remove=[Remove i1];
            end
        end
    end
end
Start=Start(setdiff(1:length(Start),Remove));
End=End(setdiff(1:length(End),Remove));
if length(Start)>length(End)
    Start=Start(1:end-1);
end

%take off artefact for entropy
Threshold=2000;
Remove=[];
for i1=1:length(Start)
    if ~isempty(find(D(D.Seizure{KeepNumber}.SelChan,Start(i1):End(i1),:)>=Threshold))
        Remove=[Remove i1];
    end
end
Start=Start(setdiff(1:length(Start),Remove));
End=End(setdiff(1:length(End),Remove));

%Local correction according to 2nd measure (ideally mean spike rate)
if CombineFlag
    for i1=1:length(Start)
        if S2(Start(i1))>=ThreshDetect2
            Start(i1)=max(find(S2(1:Start(i1))<ThreshDetect2))+1;
        else
            Start(i1)=Start(i1)+min(find(S2(Start(i1)+1:end)>=ThreshDetect2));
        end
    end
    for i1=1:length(End)
        if S2(End(i1))<ThreshDetect2
            End(i1)=max(find(S2(1:End(i1))>=ThreshDetect2));
        else
            End(i1)=End(i1)+min(find(S2(End(i1)+1:end)<ThreshDetect2))-1;
        end
    end
end

%combine Amplitude and DC
Remove=[];
if CombineAandDC
    S2=D.Seizure{NumberDC}.data;
    for i1=1:length(Start)
        if isempty(find(S2(Start(i1):End(i1))>=ThreshDetect2))
            Remove=[Remove i1];
        end
    end
    Start=Start(setdiff(1:length(Start),Remove));
    End=End(setdiff(1:length(End),Remove));
end

%use frequencies to separate seizures and somnolences
if UseFrequency && Somnolence
    Remove=[];
    Somnolence=[];
    MaxFreq=12;
    MinFreq=5.8;
    for i1=1:length(Start)
        Freq=ImaGIN_time2freq(Time(Start(i1):End(i1)));
        power=abs(fft(D(D.Seizure{KeepNumber}.SelChan,Start(i1):End(i1),:)));
        FreqPrinc=Freq(min(find(power==max(power))));
        if  FreqPrinc>=MaxFreq || FreqPrinc<=MinFreq
            Remove=[Remove i1];
            Somnolence=[Somnolence i1];
        end
    end
    SomnoStart=Start(intersect(1:length(Start),Somnolence));
    SomnoEnd=End(intersect(1:length(End),Somnolence));
    Start=Start(setdiff(1:length(Start),Remove));
    End=End(setdiff(1:length(End),Remove));
elseif UseFrequency && ~Somnolence
    Remove=[];
   MaxFreq=12;
    MinFreq=5.8;
    for i1=1:length(Start)
        Freq=ImaGIN_time2freq(Time(Start(i1):End(i1)));
        power=abs(fft(D(D.Seizure{KeepNumber}.SelChan,Start(i1):End(i1),:)));
        FreqPrinc=Freq(min(find(power==max(power))));
        if  FreqPrinc>=MaxFreq || FreqPrinc<=MinFreq
            Remove=[Remove i1];
        end
    end
    Start=Start(setdiff(1:length(Start),Remove));
    End=End(setdiff(1:length(End),Remove));
end

%add events
Index=D.Seizure{KeepNumber}.Name(end);
NewName{1}=['Seizure' Index 'Start'];
NewName{2}=['Seizure' Index 'End'];
Timing{1}=Start;
Timing{2}=End;
if Somnolence
    NewName{3}=['Somnolence' Index 'Start'];
    NewName{4}=['Somnolence' Index 'End'];
    Timing{3}=SomnoStart;
    Timing{4}=SomnoEnd;
end

Events=events(D);
if(~isempty(Events))
    try
        Events(1).type;
    catch
        Events=Events{1};
    end
end
NeventOld=size(Events,2);
for i1=1:size(Timing,2)
    for i2=1:length(Timing{i1})
        t=Time(Timing{i1});
        Events(i2+NeventOld).type=NewName(i1);
        Events(i2+NeventOld).value=i2+NeventOld;
        Events(i2+NeventOld).time=t(i2);
    end
    NeventOld=NeventOld+length(Timing{i1});
end

D=events(D,1,Events);

save(D);
spm('Pointer', 'Arrow');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_SeizureSeizureDelete(D,P,S)

NSeizure=length(D.Seizure);
if NSeizure>1
    str='';
    for i1=1:NSeizure
        str=[str num2str(i1) '|'];
    end
    str=str(1:end-1);
    KeepNumber = str2num(spm_input('Delete seizure number ',1,str));
else
    KeepNumber=1;
end

Index=D.Seizure{KeepNumber}.Name(end);
D.Seizure=D.Seizure(setdiff(1:NSeizure,KeepNumber));
Event=D.events;

%remove events
if isfield(D.events,'type')
    NewName{1}=['Seizure' Index 'Start'];
    NewName{2}=['Seizure' Index 'End'];
    for i1=1:2
        nEvent=size(Event,2);
        n=nEvent;
        i=1;
        while n>0
            if strmatch(NewName{i1},strvcat(Event(i).type))
                Event(i:nEvent-1)=Event(i+1:nEvent);
                Event(nEvent)=[];
                i=i-1;
                nEvent=nEvent-1;
            end
            n=n-1;
            i=i+1;
        end
    end
end

D=events(D,1,Event);
save(D);



