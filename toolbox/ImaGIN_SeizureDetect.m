function ImaGIN_SeizureDetect(S)
% compute instantaneous power with time-windowed fft to detect seizures
% FORMAT D = ImaGIN_SeizureDetect(S)

% try
%     Create=S.Create;
% catch
%     Flag = spm_input('Create template of seizures ',1,'Yes|No');
%     switch Flag
%         case 'Yes'
%             Create=1;
%         case 'No'
%             Create=0;
%     end
% end

% if Create
%     try
%         ImaGIN_CreateTemplate(S);
%     catch
%         ImaGIN_CreateTemplate;
%     end
% else

% try
%     ImaGIN_SeizureDetectMain(S);
% catch
%     ImaGIN_SeizureDetectMain;
% end
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ImaGIN_SeizureDetectMain(S)

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
    %     Job = spm_input('Seizure manipulation
    %     ',1,'Create|Replace|Events|Delete');
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

% if isfield(D,'time')
%     time=D.time;
% else
%     time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
%     time=time-D.timeonset;
% end
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

% try
%     Start = S.Start;
% catch
%     Start=spm_input('Start of analysis window [sec]', '+1', 'r',Time(1));
% end
% try
%     End = S.End;
% catch
%     End=spm_input('End of analysis window [sec]', '+1', 'r',Time(end));
% end

try
    Method = S.Method;
catch
%     Method = spm_input('Measure for seizure detection ',1,'Spike|AR|PowerAR|Power|Power2|Entropy');
    Method = spm_input('Measure for seizure detection ',1,'SVD|Spike|DC|Entropy|Amplitude');
end

switch Method
    case 'AR'

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
       
        
        
%         try
%             Mfactor=S.Mfactor;
%         catch
%             Mfactor = spm_input('Number of oscillations ', '+1', 'i',4,1);
%         end

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

    case 'Power'

        try
            FreqBand=S.frequencies;
        catch
            FreqBand = spm_input('Frequency band (Hz)', '+1', 'r', '', [1, 2]);
        end
% 
%         try
%             Freq2=S.Freq2;
%         catch
%             Flag = spm_input('Normalise to other frequencies ','+1','Yes|No');
%             switch Flag
%                 case 'Yes'
%                     Freq2 = spm_input('Frequencies of non interest (Hz)', '+1', 'r', '', [1, inf]);
%                 case 'No'
%                     Freq2=[];
%             end
%         end

        try
            Baseline = S.Baseline;
        catch
            Baseline = spm_select(inf, '\.mat$', 'Select EEG baseline mat file');
        end
        if isempty(Baseline)
            Dbase=[];
        else
            if size(Baseline,1)==1
                Dbase=spm_eeg_load(Baseline);
            else
                for i1=1:size(Baseline,1)
                    Dbase{i1}=spm_eeg_load(deblank(Baseline(i1,:)));
                end
            end
        end

        try
            Template = S.Template;
        catch
            Template = spm_select(inf, '\.mat$', 'Select EEG seizure mat file');
        end
        if size(Template,1)==1
            Dref=spm_eeg_load(Template);
        else
            for i1=1:size(Template,1)
                Dref{i1}=spm_eeg_load(deblank(Template(i1,:)));
            end
        end

        
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',1);
        end
        
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
        
        
%         try
%             Mfactor=S.Mfactor;
%         catch
%             Mfactor = spm_input('Number of oscillations ', '+1', 'i',4,1);
%         end

        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
%         D.Seizure{KeepNumber}.Mfactor=Mfactor;
        D.Seizure{KeepNumber}.FreqBand=FreqBand;

    case 'PowerAR'


        try
            Baseline = S.Baseline;
        catch
            Baseline = spm_select(inf, '\.mat$', 'Select EEG baseline mat file');
        end
        if isempty(Baseline)
            Dbase=[];
        else
            if size(Baseline,1)==1
                Dbase=spm_eeg_load(Baseline);
            else
                for i1=1:size(Baseline,1)
                    Dbase{i1}=spm_eeg_load(deblank(Baseline(i1,:)));
                end
            end
        end

        try
            Template = S.Template;
        catch
            Template = spm_select(inf, '\.mat$', 'Select EEG seizure mat file');
        end
        if size(Template,1)==1
            Dref=spm_eeg_load(Template);
        else
            for i1=1:size(Template,1)
                Dref{i1}=spm_eeg_load(deblank(Template(i1,:)));
            end
        end

        
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2);
        end
        
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
        
        
%         try
%             Mfactor=S.Mfactor;
%         catch
%             Mfactor = spm_input('Number of oscillations ', '+1', 'i',4,1);
%         end

        %Save parameter analysis
        
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
%         D.Seizure{KeepNumber}.Mfactor=Mfactor;

    case 'Power2'

        try
            FreqBand=S.frequencies;
        catch
            FreqBand = spm_input('Frequency band (Hz)', '+1', 'r', '', [1, 2]);
        end
        
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        end
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',1);
        end
        
        Start=min(TimeWindow)-TimeWindowWidth;
        End=max(TimeWindow)+TimeWindowWidth;
        
        
%         try
%             Mfactor=S.Mfactor;
%         catch
%             Mfactor = spm_input('Number of oscillations ', '+1', 'i',4,1);
%         end

        %Save parameter analysis
        
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.TimeWindowWidth=TimeWindowWidth;
%         D.Seizure{KeepNumber}.Mfactor=Mfactor;
        D.Seizure{KeepNumber}.FreqBand=FreqBand;

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

% try
%     Bin=S.Bin;
% catch
%     Flag = spm_input('Detect start and end of seizures ','+1','Yes|No');
%     switch Flag
%         case 'Yes'
%             Bin=1;
%         case 'No'
%             Bin=0;
%     end
% end
Bin=0;

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
        ThreshDetect = spm_input('Threshold detection (>0, 0=auto) ', '+1', 'r',0,1);
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
        
        
%         for i2=1:size(Data,1)
%             Data(i2,:)=ImaGIN_bandpassFilter(Data(i2,:),fsample(D),7,90);
%         end

        
%         [a,b]=butter(10,7*2/fsample(D));
%         b=fir1(4,7*2/fsample(D),'high');
        
        
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            
            data=Data(:,win).*(ones(size(Data,1),1)*hanning(length(win))');
            
            data=ImaGIN_Normalisation(data,2,[]);
%             data=filter(b,1,data,[],2);
%             for i2=1:size(data,1)
%                 data(i2,:)=conv(data(i2,:),b,'same');
%             end
            for i2=1:size(data,1)
                data(i2,:)=ImaGIN_bandpassFilter(data(i2,:),fsample(D),7,90);
            end
            data=ImaGIN_Normalisation(data,2,[]);
            [u,s,v]=svd(data',0);
            s=diag(s);
            
            Seizure(i1)=sum(s(1:ceil(length(s)/2)).^2)/sum(s.^2);
%             Seizure(i1)=s(1).^2/sum(s.^2);
            
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
            
%             %maxima locaux
%             DataMax= find(data>ImaGIN_shift(data,1)&data>ImaGIN_shift(data,-1));
%             %minima locaux
%             DataMin= find(data<ImaGIN_shift(data,1)&data<ImaGIN_shift(data,-1));
            
                
            
%             dmin=sort(data(DataMin));
%             dmax=sort(data(DataMax),'descend');
            dmin2=sort(data);
            dmax2=sort(data,'descend');

            data=sort(data);
            data=data(round(0.3*length(data)):round(0.7*length(data)));

            if PeakPos
%                 Seizure(i1)=(mean(dmax(1:ceil(0.1*length(dmax))))+mean(dmin(1:ceil(0.1*length(dmin)))))/std(data);
%                 Seizure(i1)=(mean(dmax(find(dmax>-dmin(ceil(0.1*length(dmin))))))+mean(dmin(1:ceil(0.1*length(dmin)))))/std(data);
%                 Seizure(i1)=(sum(dmax2(find(dmax2>-dmin(ceil(0.1*length(dmin))))))+sum(dmin2(find(dmin2<min(ceil(0.1*length(dmin)))))))/std(data);
                Seizure(i1)=(sum(dmax2(find(dmax2>-dmin2(ceil(0.1*length(dmin2))))))+sum(dmin2(find(dmin2<dmin2(ceil(0.1*length(dmin2)))))))/std(data);
            else
%                Seizure(i1)=-(mean(dmax(1:ceil(0.1*length(dmax))))+mean(dmin(1:ceil(0.1*length(dmin)))))/std(data);
%                 Seizure(i1)=-(mean(dmax(1:ceil(0.1*length(dmax))))+mean(dmin(find(dmin<-dmax(ceil(0.1*length(dmax)))))))/std(data);
%                 Seizure(i1)=-(sum(dmax2(find(dmax2>-dmin(ceil(0.1*length(dmin))))))+sum(dmin2(find(dmin2<min(ceil(0.1*length(dmin)))))))/std(data);
                Seizure(i1)=-(sum(dmax2(find(dmax2>dmax2(ceil(0.1*length(dmax2))))))+sum(dmin2(find(dmin2<-dmax2(ceil(0.1*length(dmax2)))))))/std(data);
            end

                        
            
            
%             Th=min(abs([dmin(ceil(0.1*length(dmin))) dmax(ceil(0.1*length(dmax)))]));
%             Seizure(i1)=abs(abs(sum(dmin<-Th))-abs(sum(dmax>Th)))/std(Data(win));
           
            
            
%             dmin=sort(Data(win(DataMin)));
%             dmax=sort(Data(win(DataMax)));
%             Seizure(i1)=mean([dmin(ceil(length(dmin)/2):ceil(length(dmin)/2)) dmax(ceil(length(dmax)/2):ceil(length(dmax)/2))])/std(Data(win));
        
%             Seizure(i1)=(median(Data(win))-mean(Data(win)))/std(Data(win));
%             Seizure(i1)=median(Data(win))/std(Data(win));
        end
%         if mean(Seizure)<0
%             Seizure=-Seizure;
%         end

        tmp=sort(abs(Seizure(find(Seizure<0))));
        tmp=tmp(round(0.95*length(tmp)));
        Seizure=abs(Seizure)./tmp;


        S = interp1(TimeWindow,Seizure,Time,'linear');
    
    
    case 'Amplitude'        %search for spike amplitude

        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);

        
%         %Apply moving average to get smooth time series
%         windowSize = round(D.fsample/30);  %33ms time window
%         Data=filter(ones(1,windowSize)/windowSize,1,Data);
        if isempty(Baseline)
            Data=ImaGIN_Normalisation(Data,2,[]);
        else
            Bsl=[];
            for i=1:length(Baseline)
                Bsl(i)=indsample(D,Baseline(i));
            end
            Data=ImaGIN_Normalisation(Data,2,Bsl(1):Bsl(2));
        end
        SpikeWidthSamples=SpikeWidth*D.fsample;
        
        %find local Maxima & Minima
        MinLoc=zeros(size(Data));
        MaxLoc=zeros(size(Data));
        for i1=3:length(Data)-2
            if Data(i1)==Data(i1-1)&Data(i1)==Data(i1+1)
                if Data(i1)-Data(i1-2)<=0&Data(i1)-Data(i1+2)<=0
                    MinLoc(i1)=1;
                end
                if Data(i1)-Data(i1-2)>=0&Data(i1)-Data(i1+2)>=0
                    MaxLoc(i1)=1;
                end
            else
                if Data(i1)-Data(i1-1)<=0&Data(i1)-Data(i1+1)<=0
                    MinLoc(i1)=1;
                end
                if Data(i1)-Data(i1-1)>=0&Data(i1)-Data(i1+1)>=0
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
            
%             %Compute differences
%             if IndexMax(1)>IndexMin(1)
%                 IndexMin2=ImaGIN_shift(IndexMin,-1);
%                 IndexMin2=IndexMin2(1:end-1);
%                 IndexMax2=IndexMax;
%             else
%                 IndexMax2=ImaGIN_shift(IndexMax,-1);
%                 IndexMax2=IndexMax2(1:end-1);
%                 IndexMin2=IndexMin;
%             end
%             n=min([length(IndexMax) length(IndexMin)]);
%             n2=min([length(IndexMax2) length(IndexMin2)]);
%             data=[Data(IndexMax(1:n))-Data(IndexMin(1:n)) Data(IndexMax2(1:n2))-Data(IndexMin2(1:n2))];
%             data2=sort(data);
%             %Normalisation constant
%             S=std(data2(1:round(1*length(data))));
            
            for i2=1:length(IndexMax)
                if Data(IndexMax(i2))>0
                    tmp=find(abs(IndexMax(i2)-IndexMin)<SpikeWidthSamples&Data(IndexMin)<0);
%                     S2=max(Data(IndexMax(i2))-Data(IndexMin(tmp)))/S;
                    S2=max(Data(IndexMax(i2))-Data(IndexMin(tmp)));
                    if S2>Seizure(i1)
                        Seizure(i1)=S2;
                    end
                end
            end
        end
                   

        S = interp1(TimeWindow,Seizure,Time,'linear');
    
    
    case 'Asymetry'        %search for signal asymetry during seizure

        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);

        
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            
            %maxima locaux
            DataMax=find(Data(win)>ImaGIN_shift(Data(win),1)&Data(win)>ImaGIN_shift(Data(win),-1));
            %minima locaux
            DataMin=find(Data(win)<ImaGIN_shift(Data(win),1)&Data(win)<ImaGIN_shift(Data(win),-1));
            dmin=sort(Data(win(DataMin)));
            dmax=sort(Data(win(DataMax)),'descend');

            Seizure(i1)=abs(abs(mean(dmin(1:ceil(0.1*length(dmin)))))-abs(mean(dmax(1:ceil(0.1*length(dmax))))))/std(Data(win));

%             dmin=sort(Data(win(DataMin)));
%             dmax=sort(Data(win(DataMax)));
%             Seizure(i1)=mean([dmin(ceil(length(dmin)/2):ceil(length(dmin)/2)) dmax(ceil(length(dmax)/2):ceil(length(dmax)/2))])/std(Data(win));
        
%             Seizure(i1)=(median(Data(win))-mean(Data(win)))/std(Data(win));
%             Seizure(i1)=median(Data(win))/std(Data(win));
        end
        if mean(Seizure)<0
            Seizure=-Seizure;
        end
        S = interp1(TimeWindow,Seizure,Time,'linear');
    
    
    case 'Entropy'        %permutation entropy

%         k=4;

%notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);


%         Coarse=ceil(D.fsample/256);
%         Coarse=ceil(D.fsample/512);
% %         Coarse=1;
        Coarse=S.Coarse;
        TimeDelay=1;
        Seizure=ones(length(TimeDelay),length(TimeWindow));
        Power=ones(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            t=Time(win);
            f = ImaGIN_Time2Freq(t);
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
                Seizure(i1)=ImaGIN_PermutationEntropy(Data(win),EmbeddingDimension,TimeDelay);
            elseif mean(power(FreqInterest))>mean(power(FreqNoInterest))
                win=win(1:Coarse:floor(TimeWindowWidth*D.fsample));
                Seizure(i1)=ImaGIN_PermutationEntropy(Data(win),EmbeddingDimension,TimeDelay);
%             for i2=1:length(TimeDelay)
%                 Seizure(i2,i1)=ImaGIN_PermutationEntropy(Data(win),EmbeddingDimension,TimeDelay(i2));
%             end
            end
            if isnan(Seizure(i1))
            end
        end

%         %Try to distinguish seizures and artefacts based on more powerful power for artefacts
%         y=ImaGIN_kMeansCluster(Seizure',2);
%         c=zeros(1,2);
%         c(1)=mean(y(find(y(:,end)==1)));
%         c(2)=mean(y(find(y(:,end)==2)));
%         s=zeros(1,2);
%         s(1)=std(y(find(y(:,end)==1)));
%         s(2)=std(y(find(y(:,end)==2)));
%         [c,order]=sort(c);
%         s=s(order);
% %         Thresh=(c(1)*s(2)+c(2)*s(1))/sum(s);
%         Thresh=(c(1)+c(2))/2;
%         TimeIndex=find(Seizure>Thresh);
%         Coarse=max([1 floor(length(TimeIndex)/2e3)]);
%         X=ImaGIN_Normalisation(Power(TimeIndex(1:Coarse:end)),2)';
% %         [W,M,R,Tlogl] = gmmbvl_em(X,4,4,0,0,0);
%         [W,M,R,Tlogl] = gmmbvl_em(X,2,2,0,0,0);
%         [W,order]=sort(W,'descend');
%         M=M(order);
%         R=R(order);
%         Thresh=(M(1)+5*R(1))*std(Power(TimeIndex(1:Coarse:end)))+mean(Power(TimeIndex(1:Coarse:end)));
% %         [M,order]=sort(M);
% %         W=W(order);
% %         R=R(order);
% %         Thresh=(M(1)+2*R(1))*std(Power(TimeIndex(1:Coarse:end)))+mean(Power(TimeIndex(1:Coarse:end)));
% %         x = linspace(min(X) - 3*max(R), max(X) + 3*max(R), 500 )';
% %         L = gmmbvl_em_gauss(x,M,R);
% %         L=L.*repmat(W',size(L,1),1);
% %         Thresh=x(min(find(L(:,2)>L(:,1)&[1:size(L,1)]'>find(L(:,1)==max(L(:,1))))))*std(Power(TimeIndex(1:Coarse:end)))+mean(Power(TimeIndex(1:Coarse:end)));
% 
%         %Apply threshold on power (correct structured motion artefact)
%         Seizure(find(Power>Thresh))=1;
        
        S = interp1(TimeWindow,Seizure,Time,'linear');
                    
        
    case 'AR'        %fft

        
        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);
        
        
        
        
        %         Freq1=Freq;
        win=find(Time>=TimeWindow(1)-TimeWindowWidth/2&Time<=TimeWindow(1)+TimeWindowWidth/2);
        win=win(1:floor(TimeWindowWidth*D.fsample));
        Time=Time(1:length(win));
        Freq=ImaGIN_Time2Freq(Time);
        Freq=Freq(2:floor(length(Freq)/2));

%         %Compute seizure template 
%         if iscell(Dref)
%             TimeRefTot=0;
%             for i1=1:length(Dref)
%                 if isfield(Dref{i1},'Time')
%                     TimeRef=Dref{i1}.Time;
%                 else
%                     TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
%                     TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
%                 end
%                 TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
%                 TimeRefTot=TimeRefTot+length(TimeWindowRef);
%             end
%             SeizureRef=zeros(TimeRefTot,1);
%             n=0;
%             for i1=1:length(Dref)
%                 DataRef=Dref{i1}.data(SelChan,:);
%                 if isfield(Dref{i1},'Time')
%                     TimeRef=Dref{i1}.Time;
%                 else
%                     TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
%                     TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
%                 end
%                 TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
%                 for i2=1:length(TimeWindowRef)
%                     n=n+1;disp(n)
%                     win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
%                     win=win(1:floor(TimeWindowWidth*Dref{i1}.fsample));
% %                     [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/4);
% %                     SeizureRef(n)=size(A,2)/size(A,1);
%                     
%                     [R, scale]   = arqr(ImaGIN_Normalisation(DataRef(win),2)', length(win)/4, 1);
%                     [sbc, fpe]   = arord(R, 1, 1,length(win)-length(win)/4, 2, length(win)/4);
%                     SeizureRef(n)=find(sbc==min(sbc));
%                 end
%             end
%         else
%             DataRef=Dref.data(SelChan,:);
%             %Remove baseline
%             DataRef=BBbase*DataRef;
%             if isfield(Dref,'Time')
%                 TimeRef=Dref.Time;
%             else
%                 TimeRef=0:1/Dref.fsample:(Dref.nsamples-1)/Dref.fsample;
%                 TimeRef=TimeRef-TimeRef(Dref.TimeZero);
%             end
%             TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
%             SeizureRef=zeros(length(TimeWindowRef),1);
%             for i2=1:length(TimeWindowRef)
%                 win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
%                 win=win(1:floor(TimeWindowWidth*Dref.fsample));
% %                 [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/2);
%                                     [R, scale]   = arqr(ImaGIN_Normalisation(DataRef(win),2)', length(win)/4, 1);
%                     [sbc, fpe]   = arord(R, 1, 1,length(win)-length(win)/4, 2, length(win)/4);
%                     SeizureRef(n)=find(sbc==min(sbc));
% 
%                 SeizureRef(n)=size(A,2)/size(A,1);
%             end
%         end
% 
%         
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(Data(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
            [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
%             Seizure(i1)=length(find(diff(find(diff(squeeze(abs(S)))>0))>1));
            AC=Autocorrelation([fliplr(squeeze(abs(S))') squeeze(abs(S))'],0.2);
%             Seizure(i1)=length(find(AC>ImaGIN_shift(AC,1)&AC>ImaGIN_shift(AC,-1)&AC<1));
%             Seizure(i1)=sum(AC(find(AC>ImaGIN_shift(AC,1)&AC>ImaGIN_shift(AC,-1)&AC<1)));
            Seizure(i1)=length(find(AC(find(AC>ImaGIN_shift(AC,1)&AC>ImaGIN_shift(AC,-1)&AC<1))>0.2));
            if mod(i1,10)==0;disp([i1 length(TimeWindow)]);plot(Seizure(1300:i1));drawnow;end
        end
        S = interp1(TimeWindow,Seizure,Time,'linear');

 
    case 'Power'        %fft

        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);

        
        
        
        %         Freq1=Freq;
        win=find(Time>=TimeWindow(1)-TimeWindowWidth/2&Time<=TimeWindow(1)+TimeWindowWidth/2);
        win=win(1:floor(TimeWindowWidth*D.fsample));
        Time=Time(1:length(win));
        Freq=ImaGIN_Time2Freq(Time);
        IndexFreqTot=2:floor(length(win)/2);
        IndexFreq=zeros(1,length(FreqBand));
        for i1=1:length(FreqBand)
            IndexFreq(i1)=min(find(abs(Freq(IndexFreqTot)-FreqBand(i1))==min(abs(Freq(IndexFreqTot)-FreqBand(i1)))));
        end
        IndexFreq=IndexFreq(1):IndexFreq(2);
        IndexFreq2=setdiff(1:length(IndexFreqTot),IndexFreq);
        
%         IndexFreq1=unique(IndexFreq1);
%         if ~isempty(Freq2)
%             IndexFreq2=zeros(1,length(Freq2));
%             for i1=1:length(Freq2)
%                 IndexFreq2(i1)=min(find(abs(Freq-Freq2(i1))==min(abs(Freq-Freq2(i1)))));
%             end
%             IndexFreq2=unique(IndexFreq2);
%         end

        %Compute baseline template 
        if ~isempty(Dbase)
            if iscell(Dbase)
                TimeRefTot=0;
                for i1=1:length(Dbase)
                    if isfield(Dbase{i1},'Time')
                        TimeRef=Dbase{i1}.Time;
                    else
                        TimeRef=0:1/Dbase{i1}.fsample:(Dbase{i1}.nsamples-1)/Dbase{i1}.fsample;
                        TimeRef=TimeRef-TimeRef(Dbase{i1}.TimeZero);
                    end
                    TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                    TimeRefTot=TimeRefTot+length(TimeWindowRef);
                end
                SeizureRef=zeros(TimeRefTot,floor(TimeWindowWidth*Dbase{i1}.fsample));
                n=0;
                for i1=1:length(Dbase)
                    DataRef=Dbase{i1}(SelChan,:);
                    if isfield(Dbase{i1},'Time')
                        TimeRef=Dbase{i1}.Time;
                    else
                        TimeRef=0:1/Dbase{i1}.fsample:(Dbase{i1}.nsamples-1)/Dbase{i1}.fsample;
                        TimeRef=TimeRef-TimeRef(Dbase{i1}.TimeZero);
                    end
                    TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                    for i2=1:length(TimeWindowRef)
                        n=n+1;
                        win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                        win=win(1:floor(TimeWindowWidth*Dbase{i1}.fsample));
                        SeizureRef(n,:)=abs(fft(DataRef(win)));
                    end
                end
            else
                DataRef=Dbase(SelChan,:);
                if isfield(Dbase,'Time')
                    TimeRef=Dbase.Time;
                else
                    TimeRef=0:1/Dbase.fsample:(Dbase.nsamples-1)/Dbase.fsample;
                    TimeRef=TimeRef-TimeRef(Dbase.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                SeizureRef=zeros(length(TimeWindowRef),floor(TimeWindowWidth*Dbase.fsample));
                for i2=1:length(TimeWindowRef)
                    win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                    win=win(1:floor(TimeWindowWidth*Dbase.fsample));
                    SeizureRef(i2,:)=abs(fft(DataRef(win)));
                end
            end
            SeizureRef=SeizureRef(:,IndexFreqTot);
            SeizureRef=SeizureRef(:,IndexFreq);
%             SeizureRefMean=mean(SeizureRef)';
%             SeizureRefStd=std(SeizureRef)';
            SeizureRefMean=0;
            SeizureRefStd=1;
%             BBbase=1;
%             %take components to explain 99% of variance
            [u,s,v]=svd(SeizureRef);
            ss=diag(s);
            ss=ss.^2;
            Variance=cumsum(ss)/sum(ss);
            NCompo=min(find(Variance>0.99));
%             NCompo=1;
            BaseRef=v(:,1:NCompo);
            InvBaseRef=pinv(BaseRef);
            BBbase=eye(length(ss))-BaseRef*InvBaseRef;
        else
            BBbase=1;
            SeizureRefMean=0;
            SeizureRefStd=1;
        end


        %Compute seizure template 
        if iscell(Dref)
            TimeRefTot=0;
            for i1=1:length(Dref)
                if isfield(Dref{i1},'Time')
                    TimeRef=Dref{i1}.Time;
                else
                    TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
                    TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                TimeRefTot=TimeRefTot+length(TimeWindowRef);
            end
            SeizureRef=zeros(TimeRefTot,length(IndexFreqTot));
            n=0;
            for i1=1:length(Dref)
                DataRef=Dref{i1}(SelChan,:);
                if isfield(Dref{i1},'Time')
                    TimeRef=Dref{i1}.Time;
                else
                    TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
                    TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                for i2=1:length(TimeWindowRef)
                    n=n+1;
                    win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                    win=win(1:floor(TimeWindowWidth*Dref{i1}.fsample));
                    tmp=abs(fft(DataRef(win)));
                    tmp=tmp(:,IndexFreqTot)';
                    %Remove baseline
%                     SeizureRef(n,:)=BBbase*tmp;
                    SeizureRef(n,:)=(BBbase*tmp-SeizureRefMean)./SeizureRefStd;
                end
            end
        else
            DataRef=Dref(SelChan,:);
            %Remove baseline
            DataRef=BBbase*DataRef;
            if isfield(Dref,'Time')
                TimeRef=Dref.Time;
            else
                TimeRef=0:1/Dref.fsample:(Dref.nsamples-1)/Dref.fsample;
                TimeRef=TimeRef-TimeRef(Dref.TimeZero);
            end
            TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
            SeizureRef=zeros(length(TimeWindowRef),length(IndexFreqTot));
            for i2=1:length(TimeWindowRef)
                win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                win=win(1:floor(TimeWindowWidth*Dref.fsample));
                tmp=abs(fft(DataRef(win)));
                tmp=tmp(:,IndexFreqTot)';
                %Remove baseline
%                 SeizureRef(i2,:)=BBbase*tmp;
                    SeizureRef(i2,:)=(BBbase*tmp-SeizureRefMean)./SeizureRefStd;
            end
        end
        
        SeizureRef=SeizureRef(:,IndexFreq);
        
%         %logarithm scale + remove 1/f component
%         SeizureRef2=log(SeizureRef);
%         for i1=1:size(SeizureRef2,1)
%             SeizureRef2(i1,:)=detrend(SeizureRef2(i1,:),1);
%         end
%         SeizureRef=SeizureRef2;
        
        %take 3 components
        [u,s,v]=svd(SeizureRef);
        ss=diag(s);
        ss=ss.^2;
        Variance=cumsum(ss)/sum(ss);
        NCompo=min(find(Variance>0.99));
%         NCompo=3;
        BaseRef=v(:,1:NCompo);
        InvBaseRef=pinv(BaseRef);
        BB=eye(length(ss))-BaseRef*InvBaseRef;
        BB2=BaseRef*InvBaseRef;

        
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            d=abs(fft(Data(win)));
            d=d(IndexFreqTot)';
            %Remove baseline
            d=BBbase*d;
            
            %Normalise
            d=(d-SeizureRefMean)./SeizureRefStd;

%             %logarithm scale + remove 1/f component
%             d=detrend(log(d),1);
            
            
%             tmp=zeros(1,size(SeizureRef,1));
%             for i2=1:size(SeizureRef,1)
%                 cc=corrcoef(d,SeizureRef(i1,:));
%                 tmp(i2)=cc(2);
%             end
%             Seizure2(i1)=max(tmp);

%             tmp=BB*d(IndexFreq);
%             Seizure(i1)=norm(d)/norm(tmp);
            tmp=BB2*d(IndexFreq);
            Seizure(i1)=norm(tmp);
            
            
%             if ~isempty(Freq2)
%                 Seizure(i1)=mean(d(IndexFreq1))/mean(d(IndexFreq2));
%             else
%                 Seizure(i1)=mean(d(IndexFreq1));
%             end
        end
        S = interp1(TimeWindow,Seizure,Time,'linear');

    
    case 'PowerAR'        %fft

        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);

        
        
        %         Freq1=Freq;
        win=find(Time>=TimeWindow(1)-TimeWindowWidth/2&Time<=TimeWindow(1)+TimeWindowWidth/2);
        win=win(1:floor(TimeWindowWidth*D.fsample));
        Time=Time(1:length(win));
        Freq=ImaGIN_Time2Freq(Time);
        Freq=Freq(2:floor(length(Freq)/2));
        

        %Compute baseline template 
        if ~isempty(Dbase)
            if iscell(Dbase)
                TimeRefTot=0;
                for i1=1:length(Dbase)
                    if isfield(Dbase{i1},'Time')
                        TimeRef=Dbase{i1}.Time;
                    else
                        TimeRef=0:1/Dbase{i1}.fsample:(Dbase{i1}.nsamples-1)/Dbase{i1}.fsample;
                        TimeRef=TimeRef-TimeRef(Dbase{i1}.TimeZero);
                    end
                    TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                    TimeRefTot=TimeRefTot+length(TimeWindowRef);
                end
                SeizureRef=zeros(TimeRefTot,floor(TimeWindowWidth*Dbase{i1}.fsample));
                n=0;
                for i1=1:length(Dbase)
                    DataRef=Dbase{i1}(SelChan,:);
                    if isfield(Dbase{i1},'Time')
                        TimeRef=Dbase{i1}.Time;
                    else
                        TimeRef=0:1/Dbase{i1}.fsample:(Dbase{i1}.nsamples-1)/Dbase{i1}.fsample;
                        TimeRef=TimeRef-TimeRef(Dbase{i1}.TimeZero);
                    end
                    TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                    for i2=1:length(TimeWindowRef)
                        n=n+1;
                        win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                        win=win(1:floor(TimeWindowWidth*Dbase{i1}.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
            [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
                        SeizureRef(n,:)=abs(squeeze(S));
                    end
                end
            else
                DataRef=Dbase(SelChan,:);
                if isfield(Dbase,'Time')
                    TimeRef=Dbase.Time;
                else
                    TimeRef=0:1/Dbase.fsample:(Dbase.nsamples-1)/Dbase.fsample;
                    TimeRef=TimeRef-TimeRef(Dbase.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                SeizureRef=zeros(length(TimeWindowRef),floor(TimeWindowWidth*Dbase.fsample));
                for i2=1:length(TimeWindowRef)
                    win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                    win=win(1:floor(TimeWindowWidth*Dbase.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
            [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
                    SeizureRef(i2,:)=abs(squeeze(S));
                end
            end
            SeizureRef=SeizureRef(:,IndexFreqTot);
            SeizureRef=SeizureRef(:,IndexFreq);
%             SeizureRefMean=mean(SeizureRef)';
%             SeizureRefStd=std(SeizureRef)';
            SeizureRefMean=0;
            SeizureRefStd=1;
%             BBbase=1;
%             %take components to explain 99% of variance
            [u,s,v]=svd(SeizureRef);
            ss=diag(s);
            ss=ss.^2;
            Variance=cumsum(ss)/sum(ss);
            NCompo=min(find(Variance>0.99));
%             NCompo=1;
            BaseRef=v(:,1:NCompo);
            InvBaseRef=pinv(BaseRef);
            BBbase=eye(length(ss))-BaseRef*InvBaseRef;
        else
            BBbase=1;
            SeizureRefMean=0;
            SeizureRefStd=1;
        end




        %Compute seizure template 
        if iscell(Dref)
            TimeRefTot=0;
            for i1=1:length(Dref)
                if isfield(Dref{i1},'Time')
                    TimeRef=Dref{i1}.Time;
                else
                    TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
                    TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                TimeRefTot=TimeRefTot+length(TimeWindowRef);
            end
            SeizureRef=zeros(TimeRefTot,length(Freq));
            n=0;
            for i1=1:length(Dref)
                DataRef=Dref{i1}(SelChan,:);
                if isfield(Dref{i1},'Time')
                    TimeRef=Dref{i1}.Time;
                else
                    TimeRef=0:1/Dref{i1}.fsample:(Dref{i1}.nsamples-1)/Dref{i1}.fsample;
                    TimeRef=TimeRef-TimeRef(Dref{i1}.TimeZero);
                end
                TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
                for i2=1:length(TimeWindowRef)
                    n=n+1;
                    win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                    win=win(1:floor(TimeWindowWidth*Dref{i1}.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
                        [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
                    tmp=abs(squeeze(S));
                    %Remove baseline
%                     SeizureRef(n,:)=BBbase*tmp;
                    SeizureRef(n,:)=(BBbase*tmp-SeizureRefMean)./SeizureRefStd;
                end
            end
        else
            DataRef=Dref(SelChan,:);
            %Remove baseline
            DataRef=BBbase*DataRef;
            if isfield(Dref,'Time')
                TimeRef=Dref.Time;
            else
                TimeRef=0:1/Dref.fsample:(Dref.nsamples-1)/Dref.fsample;
                TimeRef=TimeRef-TimeRef(Dref.TimeZero);
            end
            TimeWindowRef=TimeWindowWidth/2+TimeRef(1):0.2:TimeRef(end)-TimeWindowWidth/2;
            SeizureRef=zeros(length(TimeWindowRef),length(Freq));
            for i2=1:length(TimeWindowRef)
                win=find(TimeRef>=TimeWindowRef(i2)-TimeWindowWidth/2&TimeRef<=TimeWindowRef(i2)+TimeWindowWidth/2);
                win=win(1:floor(TimeWindowWidth*Dref.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(DataRef(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
            [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
                tmp=abs(squeeze(S));
                %Remove baseline
%                 SeizureRef(i2,:)=BBbase*tmp;
                    SeizureRef(i2,:)=(BBbase*tmp-SeizureRefMean)./SeizureRefStd;
            end
        end
        
        
%         %logarithm scale + remove 1/f component
%         SeizureRef2=log(SeizureRef);
%         for i1=1:size(SeizureRef2,1)
%             SeizureRef2(i1,:)=detrend(SeizureRef2(i1,:),1);
%         end
%         SeizureRef=SeizureRef2;
        
        %take 3 components
        [u,s,v]=svd(SeizureRef);
        ss=diag(s);
        ss=ss.^2;
        Variance=cumsum(ss)/sum(ss);
        NCompo=min(find(Variance>0.99));
%         NCompo=3;
        BaseRef=v(:,1:NCompo);
        InvBaseRef=pinv(BaseRef);
        BB=eye(length(ss))-BaseRef*InvBaseRef;
        BB2=BaseRef*InvBaseRef;

        
        Seizure=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            [w,A,C,SBC,FPE,th]=arfit(ImaGIN_Normalisation(Data(win),2)',2,length(win)/4);
            Order=size(A,2)/size(A,1);
            AR2=A;
            M = size(AR2,1);
            X.A = [eye(M),-AR2];
            X.B = eye(M);
            X.C = C;
            X.D = sqrtm(X.C);
            X.datatype = 'MVAR';
            [S,COH]=ImaGIN_mvfreqz_SCOH(X.B,X.A,X.C,Freq,D.fsample);
            d=abs(squeeze(S));
            %Remove baseline
            d=BBbase*d;
            
            %Normalise
            d=(d-SeizureRefMean)./SeizureRefStd;

%             %logarithm scale + remove 1/f component
%             d=detrend(log(d),1);
            
            
%             tmp=zeros(1,size(SeizureRef,1));
%             for i2=1:size(SeizureRef,1)
%                 cc=corrcoef(d,SeizureRef(i1,:));
%                 tmp(i2)=cc(2);
%             end
%             Seizure2(i1)=max(tmp);


%             tmp=BB*d(IndexFreq);
%             Seizure(i1)=norm(d)/norm(tmp);
            tmp=BB2*d;
            Seizure(i1)=norm(tmp);
            
            
%             if ~isempty(Freq2)
%                 Seizure(i1)=mean(d(IndexFreq1))/mean(d(IndexFreq2));
%             else
%                 Seizure(i1)=mean(d(IndexFreq1));
%             end
            if mod(i1,10)==0;disp([i1 length(TimeWindow)]);plot(Seizure(1:i1));drawnow;end
        end
        S = interp1(TimeWindow,Seizure,Time,'linear');

   
        
%     case 'Power'    %Wavelet
% 
%         %Downsample data
%         Coarse=max([1 floor(D.fsample/(3*max([Freq Freq2])))]);
%         fsample=D.fsample/Coarse;
%         if Coarse>1
%             Data=ImaGIN_lowpassFilter(Data,D.fsample,fsample/2)';
%         end
%         Data=Data(:,1:Coarse:end);
%         Time2=Time(1:Coarse:end);
%         NTime=length(Time2);
% 
%         FactMod=0;
% 
%         if ~isempty(Freq2)
%             M = spm_eeg_morlet(Mfactor, 1000/fsample, Freq2);
%             Seizure2=-Inf*ones(1,size(Data,2));
%             for i0=1:length(Freq2)
%                 if FactMod~=0
%                     FrequencyMin = Freq(i0)*(1-1/(2*FactMod));
%                     FrequencyMax = Freq(i0)*(1+1/(2*FactMod));
%                     BP=ImaGIN_bandpassFilter(Data,fsample,FrequencyMin,FrequencyMax);
%                 else
%                     BP=Data;
%                 end
%                 tmp = conv(BP, M{i0});
%                 % Time shift to remove delay
%                 tmp = tmp([1:size(BP,2)] + (length(M{i0})-1)/2);
%                 A=abs(tmp);
%                 B=A;
%                 NBin2=round(length(M{i0})/2);
%                 NBin=length(M{i0});
%                 TimeWindowIndex=NBin2:ceil(NBin2/4):length(B);
%                 for i1=1:length(TimeWindowIndex)
%                     win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
%                     B(win)=min(A(win));
%                 end
%                 TimeWindowIndex=NBin:ceil(NBin/4):length(B);
%                 for i1=1:length(TimeWindowIndex)
%                     win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
%                     Seizure2(win)=max(Seizure2(win),max(B(win))*ones(size(win)));
%                 end
%             end
%         end
% %         figure(11)
% %         plot(Seizure2)
%         Seizure=-Inf*ones(1,size(Data,2));
%         M = spm_eeg_morlet(Mfactor, 1000/fsample, Freq);
%         for i0=1:length(Freq)
%             if FactMod~=0
%                 FrequencyMin = Freq(i0)*(1-1/(2*FactMod));
%                 FrequencyMax = Freq(i0)*(1+1/(2*FactMod));
%                 BP=ImaGIN_bandpassFilter(Data,fsample,FrequencyMin,FrequencyMax);
%             else
%                 BP=Data;
%             end
%             tmp = conv(BP, M{i0});
%             % Time shift to remove delay
%             tmp = tmp([1:size(BP,2)] + (length(M{i0})-1)/2);
%             A=abs(tmp);
%             if ~isempty(Freq2)
%                 A=A./Seizure2;
%             end
%             A=ImaGIN_Normalisation(A,2);
%             B=A;
%             NBin2=round(length(M{i0})/2);
%             NBin=length(M{i0});
% %             TimeWindowIndex=NBin2:ceil(NBin2/4):length(B);
% %             for i1=1:length(TimeWindowIndex)
% %                 win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
% %                 B(win)=min(A(win));
% %             end
%             TimeWindowIndex=NBin:ceil(NBin/4):length(B);
%             for i1=1:length(TimeWindowIndex)
%                 win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
%                 Seizure(win)=max(Seizure(win),max(B(win))*ones(size(win)));
%             end
%         end
%         S=Seizure-min(Seizure)+2*eps;
%         S = interp1(Time2,S,Time,'linear');

    case 'Power2'        %fft

        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);

        
        
        %         Freq1=Freq;
        win=find(Time>=TimeWindow(1)-TimeWindowWidth/2&Time<=TimeWindow(1)+TimeWindowWidth/2);
        win=win(1:floor(TimeWindowWidth*D.fsample));
        Time=Time(1:length(win));
        Freq=ImaGIN_Time2Freq(Time);
        IndexFreq=zeros(1,length(FreqBand));
        for i1=1:length(FreqBand)
            IndexFreq(i1)=min(find(abs(Freq-FreqBand(i1))==min(abs(Freq-FreqBand(i1)))));
        end
        IndexFreq=IndexFreq(1):IndexFreq(2);        %fundamental
        IndexFreq2=2*IndexFreq(1):2*IndexFreq(2);   %first harmonic
        
        
        Seizure=zeros(1,length(TimeWindow));
        Seizure2=zeros(1,length(TimeWindow));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            d=abs(fft(Data(win)));
                        
            Seizure(i1)=mean(d(IndexFreq2))/mean(d(IndexFreq));
            Seizure2(i1)=max(d(IndexFreq2))/max(d(IndexFreq));
            
        end
        S = interp1(TimeWindow,Seizure,Time,'linear');

        
        Seizure=zeros(length(Frequency),length(TimeWindow));
        Frequency=FreqBand(1):FreqBand(2);
        for i2=1:length(Frequency)
            disp([i2 length(Frequency)])
            Data1=ImaGIN_bandpassFilter(Data,D.fsample,Frequency(i2)-1,Frequency(i2)+1);
            Data2=ImaGIN_bandpassFilter(Data,D.fsample,2*Frequency(i2)-1,2*Frequency(i2)+1);
    A1=abs(hilbert(Data1));
    A2=abs(hilbert(Data2));
        for i1=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i1)-TimeWindowWidth/2&Time<=TimeWindow(i1)+TimeWindowWidth/2);
            win=win(1:floor(TimeWindowWidth*D.fsample));
            Seizure(i2,i1)=max(A2(win))/max(A1(win));
            
        end
        plot(max(Seizure(1:i2,:)));drawnow
        end
        
        
    case 'Spike'

        
        %notch filter
Wo = 50*(Time(2)-Time(1))*2;  BW = Wo/35;
Data = ImaGIN_notch(Data, Wo, BW);
        
        
        [S,Event1,Event2]=ImaGIN_SpikeDetect(Data,Time,Freq,ThreshData,ThreshCC,0,Coarse);
        D.Seizure{KeepNumber}.Event1=Event1;
        D.Seizure{KeepNumber}.Event2=Event2;
end




if Start>1
    S=[zeros(1,Start-1) S];
end
if End<D.nsamples
    S=[S zeros(1,D.nsamples-End)];
end




% Seizure=zeros(1,length(TimeWindowIndex));
% for i0=1:length(Freq)
%     BP=ImaGIN_bandpassFilter(Data,D.fsample,Freq(i0)-FreqWidth/2,Freq(i0)+FreqWidth/2);
%     BP=ImaGIN_Normalisation(BP,2);
%     A=hilbert(BP);
%     B=abs(A);
%     for i1=1:length(TimeWindowIndex)
%         win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
%         B(win)=min(abs(A(win)));
%     end
%     for i1=1:length(TimeWindowIndex)
%         win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
%         Seizure(i1)=max([Seizure(i1) max(abs(B(win)))]);
%     end
% end
% S = interp1(Time(1:Coarse:end),Seizure,Time,'linear');
% if Start>1
%     S=[zeros(1,Start-1) S];
% end
% if End<D.nsamples
%     S=[S zeros(1,D.nsamples-End)];
% end
D.Seizure{KeepNumber}.data=S;

if Bin
    if ThreshDetect==0
        %Determine the threshold with a 4 Gaussian mixture model
        Coarse=max([1 floor(length(S)/2e3)]);
        X=ImaGIN_Normalisation(S(1:Coarse:end),2)';
        [W,M,R,Tlogl] = gmmbvl_em(X,4,4,0,0,0);
        [M,order]=sort(M);
        W=W(order);
        R=R(order);
        x = linspace(min(X) - 3*max(R), max(X) + 3*max(R), 500 )';
        L = gmmbvl_em_gauss(x,M,R);
        L=L.*repmat(W',size(L,1),1);
        Thresh=x(max(find(L(:,end-1)>L(:,end))))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
        %     Thresh=x(min(find(L(:,end)>max(L(:,end))/2)))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
    else
        %     SS=sort(S);
        %     if ThreshDetect>1
        Thresh=ThreshDetect;
        %     else
        %         Thresh=SS(round(length(S)*ThreshDetect));
        %     end
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


    %     RemoveStart=[];
    %     RemoveEnd=[];
    %     for i1=1:length(Start)
    %         try
    %             if Start(i1+1)-End(i1)<SeizureInterval&End(i1)-Start(i1)<SeizureDuration&End(i1+1)-Start(i1+1)<SeizureDuration
    %                 RemoveStart=[RemoveStart i1];
    %                 RemoveEnd=[RemoveEnd i1];
    %             end
    %         end
    %     end
    %     Start=Start(setdiff(1:length(Start),RemoveStart));
    %     End=End(setdiff(1:length(End),RemoveEnd));
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
        Events(i1+NeventOld).time=Time(Timing{i1});
    end
    Events(1+NeventOld).duration=Time(Timing{i1})-Time(Timing{i1});

    D=events(D,1,Events);
end

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
    case 'Power'

        try
            Freq=S.frequencies;
        catch
            Freq = spm_input('Frequencies of interest (Hz)', '+1', 'r', '', [1, inf]);
        end

        try
            Freq2=S.Freq2;
        catch
            Flag = spm_input('Normalise to other frequencies ','+1','Yes|No');
            switch Flag
                case 'Yes'
                    Freq2 = spm_input('Frequencies of non interest (Hz)', '+1', 'r', '', [1, inf]);
                case 'No'
                    Freq2=[];
            end
        end

        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        end
        % if isempty(TimeWindow)
        %     try
        %         TimeResolution=S.TimeResolution;
        %     catch
        %         TimeResolution = spm_input('Time sampling rate [sec]', '+1', 'r',.1,1);
        %     end
        %     if isempty(TimeResolution)
        %         Coarse=1;
        %     else
        %         Coarse=max([1 round(TimeResolution*D.fsample)]);
        %     end
        % end

        try
            Mfactor=S.Mfactor;
        catch
            Mfactor = spm_input('Number of oscillations ', '+1', 'i',4,1);
        end

        %Save parameter analysis
        D.Seizure{KeepNumber}.TimeWindow=TimeWindow;
        D.Seizure{KeepNumber}.Mfactor=Mfactor;

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
% D.Seizure{KeepNumber}.Name=['Seizure ' num2str(NameNumber)];

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
    case 'Power'

        %Downsample data
        Coarse=max([1 floor(D.fsample/(3*max([Freq Freq2])))]);
        fsample=D.fsample/Coarse;
        if Coarse>1
            Data=ImaGIN_lowpassFilter(Data,D.fsample,fsample/2);
        end
        Data=Data(:,1:Coarse:end);
        Time2=Time(1:Coarse:end);
        NTime=length(Time2);

        FactMod=0;

        if ~isempty(Freq2)
            M = spm_eeg_morlet(Mfactor, 1000/fsample, Freq2);
            Seizure2=-Inf*ones(1,size(Data,2));
            for i0=1:length(Freq2)
                if FactMod~=0
                    FrequencyMin = Freq(i0)*(1-1/(2*FactMod));
                    FrequencyMax = Freq(i0)*(1+1/(2*FactMod));
                    BP=ImaGIN_bandpassFilter(Data,fsample,FrequencyMin,FrequencyMax);
                else
                    BP=Data;
                end
                tmp = conv(BP, M{i0});
                % time shift to remove delay
                tmp = tmp([1:size(BP,2)] + (length(M{i0})-1)/2);
                A=abs(tmp);
                B=A;
                NBin2=round(length(M{i0})/2);
                NBin=length(M{i0});
                TimeWindowIndex=NBin2:ceil(NBin2/4):length(B);
                for i1=1:length(TimeWindowIndex)
                    win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
                    B(win)=min(A(win));
                end
                TimeWindowIndex=NBin:ceil(NBin/4):length(B);
                for i1=1:length(TimeWindowIndex)
                    win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
                    Seizure2(win)=max(Seizure2(win),max(B(win))*ones(size(win)));
                end
            end
        end
%         figure(11)
%         plot(Seizure2)
        Seizure=-Inf*ones(1,size(Data,2));
        M = spm_eeg_morlet(Mfactor, 1000/fsample, Freq);
        for i0=1:length(Freq)
            if FactMod~=0
                FrequencyMin = Freq(i0)*(1-1/(2*FactMod));
                FrequencyMax = Freq(i0)*(1+1/(2*FactMod));
                BP=ImaGIN_bandpassFilter(Data,fsample,FrequencyMin,FrequencyMax);
            else
                BP=Data;
            end
            tmp = conv(BP, M{i0});
            % time shift to remove delay
            tmp = tmp([1:size(BP,2)] + (length(M{i0})-1)/2);
            A=abs(tmp);
            if ~isempty(Freq2)
                A=A./Seizure2;
            end
            A = ImaGIN_Normalisation(A,2);
            B=A;
            NBin2=round(length(M{i0})/2);
            NBin=length(M{i0});
            TimeWindowIndex=NBin2:ceil(NBin2/4):length(B);
            for i1=1:length(TimeWindowIndex)
                win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
                B(win)=min(A(win));
            end
            TimeWindowIndex=NBin:ceil(NBin/4):length(B);
            for i1=1:length(TimeWindowIndex)
                win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
                Seizure(win)=max(Seizure(win),max(B(win))*ones(size(win)));
            end
        end
        S=Seizure-min(Seizure)+2*eps;
        S = interp1(Time2,S,Time,'linear');

    case 'Spike'
        S=ImaGIN_SpikeDetect(Data,Time,Freq,ThreshData,0.5,1);
end


if Start>1
    S=[zeros(1,Start-1) S];
end
if End<D.nsamples
    S=[S zeros(1,D.nsamples-End)];
end

% Seizure=zeros(1,length(TimeWindowIndex));
% for i0=1:length(Freq)
%     BP=ImaGIN_bandpassFilter(Data,D.fsample,Freq(i0)-FreqWidth/2,Freq(i0)+FreqWidth/2);
%     BP=ImaGIN_Normalisation(BP,2);
%     A=hilbert(BP);
%     B=abs(A);
%     for i1=1:length(TimeWindowIndex)
%         win=max([1 TimeWindowIndex(i1)-NBin2]):min([NTime TimeWindowIndex(i1)+NBin2]);
%         B(win)=min(abs(A(win)));
%     end
%     for i1=1:length(TimeWindowIndex)
%         win=max([1 TimeWindowIndex(i1)-NBin]):min([NTime TimeWindowIndex(i1)+NBin]);
%         Seizure(i1)=max([Seizure(i1) max(abs(B(win)))]);
%     end
% end
% S = interp1(Time(1:Coarse:end),Seizure,Time,'linear');
% if Start>1
%     S=[zeros(1,Start-1) S];
% end
% if End<D.nsamples
%     S=[S zeros(1,D.nsamples-End)];
% end
D.Seizure{KeepNumber}.data=S;

if Bin
    if ThreshDetect==0
        %Determine the threshold with a 4 Gaussian mixture model
        Coarse=max([1 floor(length(S)/2e3)]);
        X=ImaGIN_Normalisation(S(1:Coarse:end),2)';
        [W,M,R,Tlogl] = gmmbvl_em(X,4,4,0,0,0);
        [M,order]=sort(M);
        W=W(order);
        R=R(order);
        x = linspace(min(X) - 3*max(R), max(X) + 3*max(R), 500 )';
        L = gmmbvl_em_gauss(x,M,R);
        L=L.*repmat(W',size(L,1),1);
        Thresh=x(max(find(L(:,end-1)>L(:,end))))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
        %     Thresh=x(min(find(L(:,end)>max(L(:,end))/2)))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
    else
        %     SS=sort(S);
        %     if ThreshDetect>1
        Thresh=ThreshDetect;
        %     else
        %         Thresh=SS(round(length(S)*ThreshDetect));
        %     end
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


    %     RemoveStart=[];
    %     RemoveEnd=[];
    %     for i1=1:length(Start)
    %         try
    %             if Start(i1+1)-End(i1)<SeizureInterval&End(i1)-Start(i1)<SeizureDuration&End(i1+1)-Start(i1+1)<SeizureDuration
    %                 RemoveStart=[RemoveStart i1];
    %                 RemoveEnd=[RemoveEnd i1];
    %             end
    %         end
    %     end
    %     Start=Start(setdiff(1:length(Start),RemoveStart));
    %     End=End(setdiff(1:length(End),RemoveEnd));
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

% if isfield(D,'time')
%     time=D.time;
% else
%     time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
%     time=time-D.timeonset;
% end
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
        NumberA = str2num(spm_input('Amplitude on seizure number ','+1',str));
        NumberDC = str2num(spm_input('DC on seizure number ','+1',str));
        ThreshDetect2 = spm_input('Threshold for DC measure ', '+1', 'r',1,1);
    else
        CombineAandDC=0;
    end
end


S=D.Seizure{KeepNumber}.data;

if ThreshDetect==0

    SS=S(find(~isnan(S)))';
%     EmbeddingDimension=5;
%     TimeDelay=ceil(D.fsample/EmbeddingDimension);
%     NSample = length(SS)-TimeDelay*(EmbeddingDimension-1)-1;
% EmbeddingMatrix=[1:NSample]'*ones(1,EmbeddingDimension)+ones(NSample,1)*[0:TimeDelay:((EmbeddingDimension-1)*TimeDelay)];
% SS=SS(EmbeddingMatrix);

%     %use kmeans assuming 2 clusters
%     y=ImaGIN_kMeansCluster(SS,2);
%     c=zeros(1,2);
%     c(1)=mean(y(find(y(:,end)==1)));
%     c(2)=mean(y(find(y(:,end)==2)));
%     th=(max(c)-min(SS))*[0:499]/500+min(SS);
%     C=zeros(2,length(th));
%     for i1=1:length(th)
%         tmp1=diff(find(SS>=th(i1)));
%         tmp2=find(SS<th(i1));
%         if isempty(tmp2)
%             C(2,i1)=0;
%         elseif isnan(mean(find(tmp2>1)))
%             C(2,i1)=max(C(2,i1));
%         else
%             C(2,i1)=mean(find(tmp2>1));
%         end
%         if isnan(mean(find(tmp1>1)))
%             C(1,i1)=0;
%         else
%             C(1,i1)=mean(find(tmp1>1));
%         end
%     end
%     C(1)=max(C(2,:));
%     Crit=filter(ones(1,30)/30,1,diff(C(2,:)));
%     [tmp1,tmp2]=max(Crit);
%     tmp3=tmp2+min(find(diff(Crit(tmp2+1:end))>0));
%     Thresh=th(tmp3);
    

    
%     %use kmeans assuming 2 clusters
    y=ImaGIN_kMeansCluster(SS,2);
    c=zeros(1,2);
    c(1)=mean(y(find(y(:,end)==1)));
    c(2)=mean(y(find(y(:,end)==2)));
    Thresh=min(c)+std(c)
%     Thresh=(c(1)*length(find(y(:,2)==2))+c(2)*length(find(y(:,2)==1)))/length(y);
% %     Thresh=mean(c);
%     y2=ImaGIN_kMeansCluster(y(find(y(:,2)==find(c==min(c)))),2);
%     c2=zeros(1,2);
%     c2(1)=max(y2(find(y2(:,2)==1)));
%     c2(2)=max(y2(find(y2(:,2)==2)));
%     Thresh=min(c2);

%     %Determine the threshold with a 2 Gaussian mixture model
%     Coarse=max([1 floor(length(S)/2e3)]);
%     X=ImaGIN_Normalisation(S(1:Coarse:end),2)';
%     [W,M,R,Tlogl] = gmmbvl_em(X,4,4,0,0,0);
%     [M,order]=sort(M);
%     W=W(order);
%     R=R(order);
%     x = linspace(min(X) - 3*max(R), max(X) + 3*max(R), 500 )';
%     L = gmmbvl_em_gauss(x,M,R);
%     L=L.*repmat(W',size(L,1),1);
%     Thresh=x(max(find(L(:,end-1)>L(:,end))))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
%     %     Thresh=x(min(find(L(:,end)>max(L(:,end))/2)))*std(S(1:Coarse:end))+mean(S(1:Coarse:end));
else
    %     SS=sort(S);
    %     if ThreshDetect>1
    Thresh=ThreshDetect;
    %     else
    %         Thresh=SS(round(length(S)*ThreshDetect));
    %     end
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
    %Agregate small bits (<MinSize s or 2 samples) if the nearest large bits is
    %closer than 1 s
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
%Remove very small bits (<MinSize s or 2 samples) if the nearest large bits is
%further than 1 s
Remove=find(Flag);
% Th=1;
% Remove=[];
% for i1=1:length(Start)
%     if Flag(i1)
%         tmp=min(abs(Start(i1)-End(find(Flag==0)))/D.fsample);
%         if tmp>Th
%             tmp=min(abs(End(i1)-Start(find(Flag==0)))/D.fsample);
%             if tmp>Th
%                 Remove=[Remove i1+1];
%             end
%         end
%     end
% end
Start=Start(setdiff(1:length(Start),Remove));
End=End(setdiff(1:length(End),Remove));


%     RemoveStart=[];
%     RemoveEnd=[];
%     for i1=1:length(Start)
%         try
%             if Start(i1+1)-End(i1)<SeizureInterval&End(i1)-Start(i1)<SeizureDuration&End(i1+1)-Start(i1+1)<SeizureDuration
%                 RemoveStart=[RemoveStart i1];
%                 RemoveEnd=[RemoveEnd i1];
%             end
%         end
%     end
%     Start=Start(setdiff(1:length(Start),RemoveStart));
%     End=End(setdiff(1:length(End),RemoveEnd));
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
        Freq=ImaGIN_Time2Freq(Time(Start(i1):End(i1)));
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
        Freq=ImaGIN_Time2Freq(Time(Start(i1):End(i1)));
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

if isfield(D,'time')
    time=D.time;
else
    time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
    time=time-timeonset(D);
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






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_CreateTemplate(S)


try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select EEG mat file(s) for template');
end

try
    SelChan = S.Channel;
catch
    SelChan = spm_input('Select channel(s) ', '+1', 'i');
end

for i1=1:size(DD,1)
    D=spm_eeg_load(deblank(DD(i1,:)));
    tmp=D(SelChan,:)';
    Data(i1,:)=tmp(:);
end


Spectrum=mean(abs(fft(Data,[],2)));

ratio=0.1;
Horizon=round(D.nsamples*ratio);
Index=-Horizon:Horizon;
for i1=2:size(Data,1)
    AC=Crosscorrelation(Data(1,:),Data(i1,:),ratio);
    [tmp1,tmp2]=max(AC);
    Data(i1,:)=ImaGIN_shift(Data(i1,:),Index(tmp2));
end
Data=Data(:,Horizon+1:end-Horizon);
[u,s,v]=svd(Data,0);
Data=v(:,find(diag(s)>max(diag(s))/2))';

tmp=strfind(D.fname,'_');
D=clone(D,[D.fnamedat(1:tmp(2)) 'Seizure' D.fnamedat(tmp(3):end)], [size(Data,1) size(Data,2) 1]);
D(:,:,:)=Data;

D=timeonset(D,1);
D=events(D,1,[]);
for i1=1:D.nchannels
    D.chanlabels{i1}=num2str(i1);
end

save(D);

spm('Pointer', 'Arrow');

