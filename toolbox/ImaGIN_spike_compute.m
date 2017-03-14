function ImaGIN_spike_compute(S)
% compute mean firing rate
% FORMAT D = ImaGIN_spike_compute(S)

try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select EEG mat file');
end
D=spm_eeg_load(DD);
P=spm_str_manip(DD,'h');

if isfield(D,'spike')   %Recalculate seizure
    
    NElectrode=length(D.spike.name);
    try
        SelChan = S.Channel;
    catch
        SelChan = spm_input('Select channel(s) ', 1, 'i',1:NElectrode);
    end

    try
        SelNeuro = S.Neurone;
    catch
        for i1=1:length(SelChan)
            NNeurone=max(D.spike.markers{SelChan(i1)}(:,1,1,1));
            SelNeuro{i1} = spm_input(['Select neurone(s) in channel ' num2str(SelChan(i1)) ' '], '+1', 'i',1:NNeurone);
        end
    end

    Time=time(D);
    
    try
        Instantaneous=S.Instantaneous;
    catch
        Ctype = {
            'Yes',...
            'No'};
        str   = 'Instantaneous';
        Sel   = spm_input(str, '+1', 'm', Ctype);
        if Sel==1
            Instantaneous = 1;
        else
            Instantaneous = 0;
        end
    end
    
    try
        Causal=S.Causal;
    catch
        Ctype = {
            'Yes',...
            'No'};
        str   = 'Causal';
        Sel   = spm_input(str, '+1', 'm', Ctype);
        if Sel==1
            Causal = 1;
        else
            Causal = 0;
        end
    end
    
    if Instantaneous==0;
        try
            TimeWindow=S.TimeWindow;
        catch
            TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        end
    else
        TimeWindow=[];
    end

    if isempty(TimeWindow)
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
        if Instantaneous==0;
            try
                TimeResolution=S.TimeResolution;
            catch
                TimeResolution = spm_input('Time sampling rate [sec]', '+1', 'r',.1,1);
            end
        end
    else
        Start=min(TimeWindow);
        End=max(TimeWindow);
    end
    
    if Instantaneous==0;
        try
            TimeWindowWidth=S.TimeWindowWidth;
        catch
            TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r',2,1);
        end
    else
        TimeWindowWidth=0;
        TimeResolution=0;
        try
            TimeBinDuration=S.TimeBinDuration;
        catch
            TimeBinDuration = spm_input('Timebin duration [sec]', '+1', 'r',0.05,1);
        end
    end
    
    try
        TemporalSmoothing=S.TemporalSmoothing;
    catch
        TemporalSmoothing = spm_input('Temporal smoothing [sec]', '+1', 'r',0,1);
    end
        
        
    %Save parameter analysis
    for i1=1:length(SelChan)
        D.spike.Analysis{SelChan(i1)}.SelNeuro=SelNeuro;
        D.spike.Analysis{SelChan(i1)}.TimeWindow=TimeWindow;
        D.spike.Analysis{SelChan(i1)}.TimeWindowWidth=TimeWindowWidth;
        if isempty(TimeWindow)
            D.spike.Analysis{SelChan(i1)}.TimeResolution=TimeResolution;
            D.spike.Analysis{SelChan(i1)}.Start=Start;
            D.spike.Analysis{SelChan(i1)}.End=End;
        end
    end
    
    %Compute mean firing rate
    if Instantaneous==0;
        if isempty(TimeWindow)
            TimeWindow=Start:TimeResolution:End;
        end
    end
    N=0;
    for i4=1:length(SelNeuro)
        N=N+length(SelNeuro{i4});
    end
    if Instantaneous
        Data=zeros(N,length(Time));
        Data=zeros(N,length(min(Time):TimeBinDuration:max(Time)));
    else
        Data=zeros(N,length(TimeWindow));
    end
    n=0;
    for i1=1:length(SelChan)
        %Crop the data
        Index=find(D.spike.timings{SelChan(i1)}>=Start&D.spike.timings{SelChan(i1)}<=End);
        timings=D.spike.timings{SelChan(i1)}(Index);
        markers=double(D.spike.markers{SelChan(i1)}(Index,1,1,1));


        
        if Instantaneous
            %Instantaneous firing rate for each neurone
            timingsindex=unique(indsample(D,timings));
            for i3=1:length(SelNeuro{i1})
                n=n+1;
%                 events=find(markers==SelNeuro{i1}(i3));
%                 for i4=1:length(events)
%                     if length(events)>1
%                         if Causal
%                             if i4==1
%                                 Data(n,timingsindex(events(i4)))=inv(timings(events(i4+1))-timings(events(i4)));
%                             else
%                                 Data(n,timingsindex(events(i4)))=inv(timings(events(i4))-timings(events(i4-1)));
%                             end
%                         else
%                             if i4==1
%                                 Data(n,timingsindex(events(i4)))=inv(timings(events(i4+1))-timings(events(i4)));
%                             elseif i4==length(events)
%                                 Data(n,timingsindex(events(i4)))=inv(timings(events(i4))-timings(events(i4-1)));
%                             else
%                                 Data(n,timingsindex(events(i4)))=2*inv(timings(events(i4+1))-timings(events(i4-1)));
%                             end
%                         end
%                     end
%                 end
                D.Firing.Name{n}=[D.spike.name{i1} ' - ' num2str(SelNeuro{i1}(i3))];
                Data(n,:)=hist(timings,min(Time):TimeBinDuration:max(Time))/TimeBinDuration;
            end
        else
            %Mean firing rate for each neurone
            for i3=1:length(SelNeuro{i1})
                n=n+1;
                for i2=1:length(TimeWindow)
                    win=TimeWindow(i2)+[-TimeWindowWidth TimeWindowWidth]/2;
                    Data(n,i2)=length(find(markers(find(timings>=win(1)&timings<=win(2)))==SelNeuro{i1}(i3)))/TimeWindowWidth;
                end
%                 windowSize = round(TemporalSmoothing*fsample(D));
%                 Data(n,:)= filter(ones(1,windowSize)/windowSize,1,Data(n,:));
                D.Firing.Name{n}=[D.spike.name{i1} ' - ' num2str(SelNeuro{i1}(i3))];
            end
        end
    end
    if Instantaneous
        S=zeros(size(Data,1),length(Time));
        for i1=1:size(Data,1)
            index=find(Data(i1,:));
%             S(i1,:) = interp1(Time(index),Data(i1,index),Time,'linear');
            S(i1,:) = interp1(min(Time):TimeBinDuration:max(Time),Data(i1,:),Time,'linear');
%             S(i1,:) = interp1(Time(index),Data(i1,index),Time,'cubic');
        end
    else
        S=zeros(size(Data,1),length(Time));
        for i1=1:size(Data,1)
            S(i1,:) = interp1(TimeWindow,Data(i1,:),Time,'linear');
        end
    end
    if TemporalSmoothing>0
        windowSize = round(TemporalSmoothing*fsample(D));
        for i1=1:size(Data,1)
        	S(i1,:)= filter(ones(1,windowSize)/windowSize,1,S(i1,:));
        end
    end

    D.Firing.data=S;

    save(D);

    spm('Pointer', 'Arrow');


else
    error('No spike data')
end


