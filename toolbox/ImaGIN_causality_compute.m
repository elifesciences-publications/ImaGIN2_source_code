function D = ImaGIN_causality_compute(S)
% Compute diffent measure of asymetric interactions

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
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
%
% Authors: Olivier David


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG causality setup',0);

try
    DD = S.D;
catch
    DD = spm_select(inf, 'mat', 'Select EEG mat file');
end

for i1=1:size(DD,1)
%     try
%         D=ImaGIN_causality_compute_main(deblank(DD(i1,:)),S{i1});
%     catch
%         try
            [D,S]=ImaGIN_causality_compute_main(deblank(DD(i1,:)),S);
%         catch
%             [D,S]=ImaGIN_causality_compute_main(deblank(DD(i1,:)));
%         end
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D,S]=ImaGIN_causality_compute_main(D,S)
    
TimeWindow=[];
TimeWindowWidth=[];
Freq=[];

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_load(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

% Check for arguments in D
if ~isfield(D, 'ca')
	D.ca = [];
end

try
    Pre=S.Pre;
catch
    Pre=spm_input('Prefix of new file', '+1', 's');
    S.Pre=Pre;
end

try
    D.ca.Method = S.Method;
catch
    Ctype = {
        'Crosscorrelation',...
        'H2',...
            'Generalised Synchronisation',...
            'Autoregressive Models'};
    str   = 'Causality estimation';
	Sel   = spm_input(str, 2, 'm', Ctype);
	D.ca.Method = Ctype{Sel};
    S.Method=D.ca.Method;
end

% try
%     SlidingTimeWindow = S.SlidingTimeWindow;
% catch
%     tmp=spm_input('Sliding time window ','+1','Yes|No');
%     switch tmp
%         case{'No'}
%             SlidingTimeWindow=0;
%         case{'Yes'}
%             SlidingTimeWindow=1;
%     end
% end

% if SlidingTimeWindow
    try
        TimeWindow = S.TimeWindow;
    catch
        TimeWindow = spm_input('Time window positions [sec]', '+1', 'r','[]');
        S.TimeWindow=TimeWindow;
    end

    if ~isempty(TimeWindow)
        try
            TimeWindowWidth = S.Width;
        catch
            TimeWindowWidth = ...
                spm_input('Time window width [sec]', '+1', 'r', '1', 1);
            S.Width=TimeWindowWidth;
        end
    end
% end
D.ca.Time=TimeWindow;
D.ca.Width=TimeWindowWidth;

%Frequency band selection
switch D.ca.Method
    case{'Autoregressive Models'}
        try
            D.ca.Freq=S.Freq;
        catch
            D.ca.Freq =spm_input('Frequency band [Hz]', '+1', 'r', '[6:36]');
        end
        S.Freq=D.ca.Freq;
        try
            D.ca.Baseline=S.Baseline;
        catch
            D.ca.Baseline =spm_input('Baseline period [sec]', '+1', 'r', '');
        end
        S.Baseline=D.ca.Baseline;
%         if ~isempty(D.ca.Baseline)
%             S.Baseline=str2num(D.ca.Baseline);
%         end
        if length(S.Baseline)~=2&length(S.Baseline)~=0
            error('Baseline should a 2 element vector or empty')
        end 
%         try
%             D.ca.FlagNormalise=S.Normalise;
%         catch
%             tmp=spm_input('Normalise data ','+1','Yes|No');
%             if strcmp(tmp,'Yes')
%                 D.ca.FlagNormalise=1;
%             else
%                 D.ca.FlagNormalise=0;
%             end
%         end
%         S.Normalise=D.ca.FlagNormalise;
        try
            D.ca.Surro=S.Surro;
        catch
            D.ca.Surro =spm_input('Surrogate type ','+1','Shift|FT1');
        end
        S.Surro=D.ca.Surro;
        try
            D.ca.pmax=S.pmax;
        catch
            D.ca.pmax =spm_input('Max p value ', '+1', 'r', '0.01');
        end
        S.pmax=D.ca.pmax;
    case{'Generalised Synchronisation'}
        try
            D.ca.MethodGS=S.MethodGS;
        catch
            D.ca.MethodGS=spm_input('Measure of generalised synchrony ','+1','Quian|Schiff');
        end
        S.MethodGS=D.ca.MethodGS;
        try
            D.ca.TimeHorizon=S.TimeHorizon;
        catch
            D.ca.TimeHorizon =spm_input('Time horizon [s] ', '+1', 'r', '0.005');
        end
        try
            D.ca.Surro=S.Surro;
        catch
            D.ca.Surro =spm_input('Surrogate type ','+1','Shift|FT1');
        end
        S.Surro=D.ca.Surro;
        try
            D.ca.pmax=S.pmax;
        catch
            D.ca.pmax =spm_input('Max p value ', '+1', 'r', '0.01');
        end
        S.pmax=D.ca.pmax;
    case{'H2'}
        try
            D.ca.TimeHorizon=S.TimeHorizon;
        catch
            D.ca.TimeHorizon =spm_input('Time horizon [s] ', '+1', 'r', '0.005');
        end
        try
            D.ca.Surro=S.Surro;
        catch
            D.ca.Surro =spm_input('Surrogate type ','+1','Shift|FT1');
        end
        S.Surro=D.ca.Surro;
        try
            D.ca.pmax=S.pmax;
        catch
            D.ca.pmax =spm_input('Max p value ', '+1', 'r', '0.01');
        end
        S.pmax=D.ca.pmax;
end

% NB: D.ca.channels maps directly into the data. To retrieve the position of the channel,
% use D.channels.order
try
    D.ca.channels = S.channels;
catch
    D.ca.channels = ...
	 	spm_input('Select channels', '+1', 'i', sprintf('1:%d',D.nchannels));
    S.channels=D.ca.channels;
end

%D=badchannels(D,setdiff(1:D.nchannels,D.ca.channels));

spm('Pointer', 'Watch'); drawnow;


fnamedat = deblank(D.fnamedat);
D.datatype = 'float32-le';
D.Label=D.ca.Method;

switch D.ca.Method
    case{'Crosscorrelation','Generalised Synchronisation','H2'}
        %Write new channels
        D.ca.CM=ConnectivityMatrix(length(D.ca.channels));
        tmp=D.chanlabels;
        for i1=1:size(D.ca.CM,2)
            disp(i1)
%             D=chanlabels(D,i1,[deblank(tmp{D.ca.channels(D.ca.CM(1,i1))}) '/' deblank(tmp{D.ca.channels(D.ca.CM(2,i1))})]);
            D.ca.chanlabels{i1}=[deblank(tmp{D.ca.channels(D.ca.CM(1,i1))}) '/' deblank(tmp{D.ca.channels(D.ca.CM(2,i1))})];
        end
    case{'Autoregressive Models'}
        %Write new channels
        D.ca.chanlabels=D.chanlabels(D.ca.channels);
end

D1=D;
D2=D;

for k = 1 : D.ntrials
%Write new events
if isempty(TimeWindow)
    Events=[];
else
    if isfield(D,'time')
        time=D.time;
    else
        time=[0:1/D.fsample:(D.nsamples-1)/D.fsample];
        time=time+timeonset(D);
    end
    try
        Events=D.events{k};
    catch
        Events=D.events;
    end
    Events=select_events(Events, [min(TimeWindow)-TimeWindowWidth/2 max(TimeWindow)+TimeWindowWidth/2])
    %renumber events
    %     for i1=1:length(D.events.types)
    %         while D.events.types(i1)~=i1
    %             D.events.types(find(D.events.types>i1))=D.events.types(find(D.events.types>i1))-1;
    %             D.events.code(find(D.events.code>i1))=D.events.code(find(D.events.code>i1))-1;
    %         end
    %     end
end
    switch D.ca.Method
        case{'Crosscorrelation'}
            [D1,D2,D1s,D2s]=ComputeCC(D1,D2,k,TimeWindow,TimeWindowWidth,fnamedat,Pre);
        case{'H2'}
            [D1,D1s]=ComputeH2(D1,k,TimeWindow,TimeWindowWidth,fnamedat,Pre);
        case{'Generalised Synchronisation'}
            [D1,D1s]=ComputeGS(D1,k,TimeWindow,TimeWindowWidth,fnamedat,Pre);
        case{'Autoregressive Models'}
            [D1,D1s]=ComputeDTF(D1,k,TimeWindow,TimeWindowWidth,S.Baseline,fnamedat,Pre);
    end
     D1s=events(D1s,k,Events);
     try
         D2s=events(D2s,k,Events);
     end
end     

save(D1s);
try
    D2s;
    save(D2s);
end
spm('Pointer', 'Arrow');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D2,D1s,D2s]=ComputeCC(D1,D2,k,TimeWindow,TimeWindowWidth,fnamedat,Pre)

%D1=strength of coupling
%D2=delay in sec

if isempty(TimeWindow)
    d1 = zeros(size(D1.ca.CM,2),1);
    d2 = zeros(size(D2.ca.CM,2),1);
    for i1=1:size(D1.ca.CM,2)
        tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),:,k));
        tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),:,k));
        AC=Crosscorrelation(tmp1,tmp2);
        [tmp1,tmp2]=max(abs(AC));
        tmp2=tmp2-ceil(length(AC)/2);
        d1(i1)=tmp1;
        d2(i1)=tmp2/D1.fsample;
    end
else
    d1 = zeros(size(D1.ca.CM,2),length(TimeWindow));
    d2 = zeros(size(D1.ca.CM,2),length(TimeWindow));
    if isfield(D1,'time')
        time=D1.time;
    else
        time=[0:1/D1.fsample:(D1.nsamples-1)/D1.fsample];
        time=time+timeonset(D1);
    end
    for i0=1:length(TimeWindow)
        win=find(time>=TimeWindow(i0)-TimeWindowWidth/2&time<=TimeWindow(i0)+TimeWindowWidth/2);
        for i1=1:size(D1.ca.CM,2)
            tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),win,k));
            tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),win,k));
            AC=Crosscorrelation(tmp1,tmp2);
            [tmp1,tmp2]=max(abs(AC));
            tmp2=tmp2-ceil(length(AC)/2);
            d1(i1,i0)=tmp1;
            d2(i1,i0)=tmp2/D1.fsample;
        end
    end
    D1.time=TimeWindow;
    D2.time=TimeWindow;
    [tmp TimeZero]=min(TimeWindow);
    D1=timeonset(D1,tmp);
    D2=timeonset(D2,tmp);
end

if isempty(TimeWindow)
    nsamples=1;
else
    nsamples=length(TimeWindow);
end

D1s=clone(D1,['cc' Pre '_ca1_' fnamedat],[size(D1.ca.CM,2) nsamples D1.ntrials]);
D1s(:,:,k)=d1;
D2s=clone(D2,['cc' Pre '_ca2_' fnamedat],[size(D2.ca.CM,2) nsamples D2.ntrials]);
D2s(:,:,k)=d2;

spm('Pointer', 'Arrow'); drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D1s]=ComputeGS(D1,k,TimeWindow,TimeWindowWidth,fnamedat,Pre)

%D1=strength of coupling
%D2=delay in sec

% if strcmp(D1.ca.MethodGS,'Schiff')
    TimeHorizon=ceil(D1.ca.TimeHorizon*D1.fsample);
% else
%     TimeHorizon=0;
% end


if isempty(TimeWindow)
    d1 = zeros(size(D1.ca.CM,2),1);
    D1.ca.GSX = zeros(size(D1.ca.CM,2),1);
    D1.ca.GSY = zeros(size(D1.ca.CM,2),1);
    for i1=1:size(D1.ca.CM,2)
        tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),:,k));
        tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),:,k));
        TimeDelay=NaN;
        n=1;
        while isnan(TimeDelay)
            AC=Autocorrelation(tmp1,n*0.1);
            tmp=ImaGIN_shift(AC,find(AC==1));
            TimeDelay1=min(find(tmp<0));
            AC=Autocorrelation(tmp2,n*0.1);
            tmp=ImaGIN_shift(AC,find(AC==1));
            TimeDelay2=min(find(tmp<0));
            TimeDelay=ceil(mean([TimeDelay1 TimeDelay2]));
            n=n+1;
        end
%         [ED,FNNP] = EmbeddingDimension(tmp1,TimeDelay,20);
        EmbeddingDimension=3;
%         TimeHorizon=ceil(TimeDelay/2);
        [PredictionX,PredictionY] = ImaGIN_NonlinearInterdependence(tmp1,tmp2,TimeDelay,EmbeddingDimension,D1.ca.MethodGS,TimeHorizon);
        d1(i1)=PredictionX-PredictionY;
        D1.ca.GSX(i1)=PredictionX;
        D1.ca.GSY(i1)=PredictionY;
    end
%     D1.time=[];
%     D1=rmfield(D1,'timeonset');
else
    d1 = zeros(size(D1.ca.CM,2),length(TimeWindow));
    D1.ca.GSX = zeros(size(D1.ca.CM,2),length(TimeWindow));
    D1.ca.GSY = zeros(size(D1.ca.CM,2),length(TimeWindow));
    if isfield(D1,'time')
        time=D1.time;
    else
        time=[0:1/D1.fsample:(D1.nsample-1)/D1.fsample];
        time=time+timeonset(D1);
    end
    for i0=1:length(TimeWindow)
        disp([i0 length(TimeWindow)])
        win=find(time>=TimeWindow(i0)-TimeWindowWidth/2&time<=TimeWindow(i0)+TimeWindowWidth/2);
        for i1=1:size(D1.ca.CM,2)
            tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),win,k));
            tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),win,k));
            TimeDelay=NaN;
            n=1;
            while isnan(TimeDelay)
                AC=Autocorrelation(tmp1,n*0.1);
                tmp=ImaGIN_shift(AC,find(AC==1));
                TimeDelay1=min(find(tmp<0));
                AC=Autocorrelation(tmp2,n*0.1);
                tmp=ImaGIN_shift(AC,find(AC==1));
                TimeDelay2=min(find(tmp<0));
                TimeDelay=ceil(mean([TimeDelay1 TimeDelay2]));
                n=n+1;
            end
            %         [ED,FNNP] = EmbeddingDimension(tmp1,TimeDelay,20);
            EmbeddingDimension=3;
%             TimeHorizon=ceil(TimeDelay/2);
            [PredictionX,PredictionY] = ImaGIN_NonlinearInterdependence(tmp1,tmp2,TimeDelay,EmbeddingDimension,D1.ca.MethodGS,TimeHorizon);
            D1.ca.GSX(i1,i0)=PredictionX;
            D1.ca.GSY(i1,i0)=PredictionY;
            d1(i1,i0)=PredictionX-PredictionY;
        end
    end
    D1.time=TimeWindow;
    [tmp D1.TimeZero]=min(abs(TimeWindow));
end



D1.ca.p=Inf;
n=0;
D1.ca.GSX_S = [];
D1.ca.GSY_S = [];
D1.ca.GS_S = [];
NTimeWindow=max([1 length(TimeWindow)]);
data=squeeze(D1(D1.ca.channels,:,k));
while min(D1.ca.p)>D1.ca.pmax
    n=n+1;
    D1.ca.p=[1:n*NTimeWindow]/(1+n*NTimeWindow);
    %Shift or FT1 surrogates
    switch D1.ca.Surro
        case 'Shift'
            tmp=round(1.9*D1.nsamples*(rand(1,size(data,1))-0.5));
            Surro=data;
            for i0=1:size(data,1)
                Surro(i0,:)=ImaGIN_shift(data(i0,:),tmp(i0),2);
            end
        case 'FT1'
            Surro=data;
            for i0=1:size(data,1)
                Surro(i0,:)=ImaGIN_FT1Surrogate(data(i0,:),1);
            end
    end

    %Surrogates
    GS_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    GSX_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    GSY_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    if isempty(TimeWindow)
        for i1=1:size(D1.ca.CM,2)
            tmp1=Surro(D1.ca.CM(1,i1),:);
            tmp2=Surro(D1.ca.CM(2,i1),:);
            TimeDelay=NaN;
            n=1;
            while isnan(TimeDelay)
                AC=Autocorrelation(tmp1,n*0.1);
                tmp=ImaGIN_shift(AC,find(AC==1));
                TimeDelay1=min(find(tmp<0));
                AC=Autocorrelation(tmp2,n*0.1);
                tmp=ImaGIN_shift(AC,find(AC==1));
                TimeDelay2=min(find(tmp<0));
                TimeDelay=ceil(mean([TimeDelay1 TimeDelay2]));
                n=n+1;
            end
            %         [ED,FNNP] = EmbeddingDimension(tmp1,TimeDelay,20);
            EmbeddingDimension=3;
            TimeHorizon=TimeDelay;
            [PredictionX,PredictionY] = ImaGIN_NonlinearInterdependence(tmp1,tmp2,TimeDelay,EmbeddingDimension,D1.ca.MethodGS,TimeHorizon);
            GSX_S(i1)=PredictionX;
            GSY_S(i1)=PredictionY;
            GS_S(i1)=PredictionX-PredictionY;
        end
    else
        for i0=1:length(TimeWindow)
            win=find(time>=TimeWindow(i0)-TimeWindowWidth/2&time<=TimeWindow(i0)+TimeWindowWidth/2);
            for i1=1:size(D1.ca.CM,2)
                tmp1=Surro(D1.ca.CM(1,i1),win);
                tmp2=Surro(D1.ca.CM(2,i1),win);
                TimeDelay=NaN;
                n=1;
                while isnan(TimeDelay)
                    AC=Autocorrelation(tmp1,n*0.1);
                    tmp=ImaGIN_shift(AC,find(AC==1));
                    TimeDelay1=min(find(tmp<0));
                    AC=Autocorrelation(tmp2,n*0.1);
                    tmp=ImaGIN_shift(AC,find(AC==1));
                    TimeDelay2=min(find(tmp<0));
                    TimeDelay=ceil(mean([TimeDelay1 TimeDelay2]));
                    n=n+1;
                end
                %         [ED,FNNP] = EmbeddingDimension(tmp1,TimeDelay,20);
                EmbeddingDimension=3;
                TimeHorizon=TimeDelay;
                [PredictionX,PredictionY] = ImaGIN_NonlinearInterdependence(tmp1,tmp2,TimeDelay,EmbeddingDimension,D1.ca.MethodGS,TimeHorizon);
                GSX_S(i1,i0)=PredictionX;
                GSY_S(i1,i0)=PredictionY;
                GS_S(i1,i0)=PredictionX-PredictionY;
            end
        end
    end
    D1.ca.GSX_S=cat(2,D1.ca.GSX_S,GSX_S);
    D1.ca.GSY_S=cat(2,D1.ca.GSY_S,GSY_S);
    D1.ca.GS_S=cat(2,D1.ca.GS_S,GS_S);
end




% D.scale(:, 1, 1, k) = (max(abs(reshape(d, [D.Nfrequencies*D.nsample length(D.tf.channels)])))./(2^31-1))';
% d = int32(d./repmat(D.scale(:, 1, 1, k), [1 D.Nfrequencies D.nsample]));
% 
% D2.scale(:, 1, 1, k) = (max(abs(reshape(d2, [D2.Nfrequencies*D2.nsample length(D2.tf.channels)])))./32767)';
% d2 = int16(d2./repmat(D2.scale(:, 1, 1, k), [1 D2.Nfrequencies D2.nsample]));

% Dint.scale(:, 1, 1, k) = (max(abs(reshape(dint, [Dint.Nfrequencies*Dint.nsample length(Dint.tf.channels)])))./(2^31-1))';
if isempty(TimeWindow)
    nsamples=1;
else
    nsamples=length(TimeWindow);
end
D1s=clone(D1,['gs' Pre '_ca1_' fnamedat],[D1.nchannels nsamples D1.ntrials]);
D1s(:,:,k)=d1
clear D2


% fwrite(fpd1, d1, D1.datatype);
% fwrite(fpd2, d2, D2.datatype);
% 

spm('Pointer', 'Arrow'); drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D1s]=ComputeH2(D1,k,TimeWindow,TimeWindowWidth,fnamedat,Pre)

%D1=strength of coupling
%D2=delay in sec

% if strcmp(D1.ca.MethodGS,'Schiff')
    TimeHorizon=ceil(D1.ca.TimeHorizon*D1.fsample);
% else
%     TimeHorizon=0;
% end


if isempty(TimeWindow)
    d1 = zeros(size(D1.ca.CM,2),1);
    D1.ca.Delay = zeros(size(D1.ca.CM,2),1);
    D1.ca.H2X = zeros(size(D1.ca.CM,2),1);
    D1.ca.H2Y = zeros(size(D1.ca.CM,2),1);
    for i1=1:size(D1.ca.CM,2)
        tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),:,k));
        tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),:,k));
        H2 = ImaGIN_H2Compute(tmp1,tmp2,TimeHorizon);
        tmp=abs(H2(1,:)-H2(2,:));
        index=min(find(tmp==max(tmp)));
        D1.ca.Delay(i1)=(index-(length(tmp)-1)/2)/D1.fsample;
        D1.ca.H2X(i1)=H2(1,index);
        D1.ca.H2Y(i1)=H2(2,index);
        d1(i1)=H2(1,index)-H2(2,index);
    end
else
    d1 = zeros(size(D1.ca.CM,2),length(TimeWindow));
    D1.ca.Delay = zeros(size(D1.ca.CM,2),length(TimeWindow));
    D1.ca.H2X = zeros(size(D1.ca.CM,2),length(TimeWindow));
    D1.ca.H2Y = zeros(size(D1.ca.CM,2),length(TimeWindow));
    Time=time(D1);
    for i0=1:length(TimeWindow)
        disp([i0 length(TimeWindow)])
        win=find(Time>=TimeWindow(i0)-TimeWindowWidth/2&Time<=TimeWindow(i0)+TimeWindowWidth/2);
        for i1=1:size(D1.ca.CM,2)
            tmp1=squeeze(D1(D1.ca.channels(D1.ca.CM(1,i1)),win,k));
            tmp2=squeeze(D1(D1.ca.channels(D1.ca.CM(2,i1)),win,k));
            H2 = ImaGIN_H2Compute(tmp1,tmp2,TimeHorizon);
            tmp=abs(H2(1,:)-H2(2,:));
            index=min(find(tmp==max(tmp)));
            D1.ca.Delay(i1,i0)=(index-(length(tmp)-1)/2)/D1.fsample;
            D1.ca.H2X(i1,i0)=H2(1,index);
            D1.ca.H2Y(i1,i0)=H2(2,index);
            d1(i1,i0)=H2(1,index)-H2(2,index);
        end
    end
    [tmp t]=min(TimeWindow);
    D1=timeonset(D1,tmp);
    if length(TimeWindow)>1
        D1=fsample(D1,1/(TimeWindow(2)-TimeWindow(1)));
    end
end



D1.ca.p=Inf;
n=0;
D1.ca.H2X_S = [];
D1.ca.H2Y_S = [];
D1.ca.H2_S = [];
NTimeWindow=max([1 length(TimeWindow)]);
data=squeeze(D1(D1.ca.channels,:,k));
while min(D1.ca.p)>D1.ca.pmax
    n=n+1;
    D1.ca.p=[1:n*NTimeWindow]/(1+n*NTimeWindow);
    %Shift or FT1 surrogates
    switch D1.ca.Surro
        case 'Shift'
            tmp=round(1.9*D1.nsamples*(rand(1,size(data,1))-0.5));
            Surro=data;
            for i0=1:size(data,1)
                Surro(i0,:)=ImaGIN_shift(data(i0,:),tmp(i0),2);
            end
        case 'FT1'
            Surro=data;
            for i0=1:size(data,1)
                Surro(i0,:)=ImaGIN_FT1Surrogate(data(i0,:),1);
            end
    end

    %Surrogates
    H2_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    H2X_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    H2Y_S = zeros(size(D1.ca.CM,2),NTimeWindow);
    if isempty(TimeWindow)
        for i1=1:size(D1.ca.CM,2)
            tmp1=Surro(D1.ca.CM(1,i1),:);
            tmp2=Surro(D1.ca.CM(2,i1),:);
            H2 = ImaGIN_H2Compute(tmp1,tmp2,TimeHorizon);
            tmp=abs(H2(1,:)-H2(2,:));
            index=min(find(tmp==max(tmp)));
            H2X_S(i1)=H2(1,index);
            H2Y_S(i1)=H2(2,index);
            H2_S(i1)=H2(1,index)-H2(2,index);
        end
    else
        for i0=1:length(TimeWindow)
            win=find(Time>=TimeWindow(i0)-TimeWindowWidth/2&Time<=TimeWindow(i0)+TimeWindowWidth/2);
            for i1=1:size(D1.ca.CM,2)
                tmp1=Surro(D1.ca.CM(1,i1),win);
                tmp2=Surro(D1.ca.CM(2,i1),win);
                H2 = ImaGIN_H2Compute(tmp1,tmp2,TimeHorizon);
                tmp=abs(H2(1,:)-H2(2,:));
                index=min(find(tmp==max(tmp)));
                H2X_S(i1,i0)=H2(1,index);
                H2Y_S(i1,i0)=H2(2,index);
                H2_S(i1,i0)=H2(1,index)-H2(2,index);
            end
        end
    end
    D1.ca.H2X_S=cat(2,D1.ca.H2X_S,H2X_S);
    D1.ca.H2Y_S=cat(2,D1.ca.H2Y_S,H2Y_S);
    D1.ca.H2_S=cat(2,D1.ca.H2_S,H2_S);
end




% D.scale(:, 1, 1, k) = (max(abs(reshape(d, [D.Nfrequencies*D.nsample length(D.tf.channels)])))./(2^31-1))';
% d = int32(d./repmat(D.scale(:, 1, 1, k), [1 D.Nfrequencies D.nsample]));
% 
% D2.scale(:, 1, 1, k) = (max(abs(reshape(d2, [D2.Nfrequencies*D2.nsample length(D2.tf.channels)])))./32767)';
% d2 = int16(d2./repmat(D2.scale(:, 1, 1, k), [1 D2.Nfrequencies D2.nsample]));

% Dint.scale(:, 1, 1, k) = (max(abs(reshape(dint, [Dint.Nfrequencies*Dint.nsample length(Dint.tf.channels)])))./(2^31-1))';
D1s=clone(D1,['h2' Pre '_ca1_' fnamedat], [size(d1,1) size(d1,2) D1.ntrials]);
for i1=1:length(D1.ca.chanlabels)
    D1s=chanlabels(D1s,i1,D1.ca.chanlabels{i1});
end
D1s(:,:,k)=d1;
save(D1s);
clear D2
% fwrite(fpd1, d1, D1.datatype);
% fwrite(fpd2, d2, D2.datatype);
% 

spm('Pointer', 'Arrow'); drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D1,D1s]=ComputeDTF(D1,k,TimeWindow,TimeWindowWidth,Baseline,fnamedat,Pre);

if isfield(D1,'time')
    time=D1.time;
else
    time=[0:1/D1.fsample:(D1.nsamples-1)/D1.fsample];
    time=time+timeonset(D1);
end

% if D1.ca.FlagNormalise
%     data=ImaGIN_Normalisation(squeeze(D1.data(D1.ca.channels,:,k))',1);
% else
%     data=squeeze(D1.data(D1.ca.channels,:,k))';
% end

%Amplitude normalisation according to the baseline
Baseline=find(time>=Baseline(1)&time<=Baseline(2));
if isempty(Baseline)
    data=ImaGIN_Normalisation(squeeze(D1(D1.ca.channels,:,k))',1);
%     data=squeeze(D1.data(D1.ca.channels,:,k))';
    %Find error covariance during the baseline (used for GPDC)
    [w,A,C,SBC,FPE,th]=arfit(data,2,ceil(D1.fsample/8));
%     [w,A,C,SBC,FPE,th]=arfit(data,2,ceil(D1.fsample/4),'fpe');
%     [AR2,RC,PE] = mvar(data,size(A,2)/size(A,1));
%     ErrorCov = PE(:,size(AR2,1)*size(A,2)/size(A,1)+1:size(AR2,1)*(size(A,2)/size(A,1)+1));
    ErrorCov=C;
else
    data=ImaGIN_Normalisation(squeeze(D1(D1.ca.channels,:,k))',1,Baseline);
%     data=squeeze(D1.data(D1.ca.channels,:,k))';
    %Find error covariance during the baseline (used for GPDC)
    [w,A,C,SBC,FPE,th]=arfit(data(Baseline,:),2,ceil(D1.fsample/8));
%     [w,A,C,SBC,FPE,th]=arfit(data(Baseline,:),2,ceil(D1.fsample/4),'fpe');
        
%     [AR2,RC,PE] = mvar(data(Baseline,:),size(A,2)/size(A,1));
%     ErrorCov = PE(:,size(AR2,1)*size(A,2)/size(A,1)+1:size(AR2,1)*(size(A,2)/size(A,1)+1));
    ErrorCov=C;
end

%Bandpass filter
% data=ImaGIN_bandpassFilter(data,D1.fsample,D1.ca.Freq(1),D1.ca.Freq(end));
%Bandpass filter (Butterworth)
% ord=buttord([0.8*D1.ca.Freq(1)]/(0.5*D1.fsample),[D1.ca.Freq(1)]/(0.5*D1.fsample),3,40);
% [b,a] = butter(10,[D1.ca.Freq(1)]/(0.5*D1.fsample),'high');
% data=filter(b,a,data);
% [b,a] = butter(10,[D1.ca.Freq(end)]/(0.5*D1.fsample),'low');
% data=filter(b,a,data);

if isempty(TimeWindow)
    %Model order
    [w,A,C,SBC,FPE,th]=arfit(data,2,ceil(D1.fsample/8));
%     [w,A,C,SBC,FPE,th]=arfit(data,2,ceil(D1.fsample/4),'fpe');
    D1.ca.Order=size(A,2)/size(A,1);
%     [AR2,RC,PE] = mvar(data,D1.ca.Order);
    AR2=A;
    M = size(AR2,1);
    X.A = [eye(M),-AR2];
    X.B = eye(M);
%     X.C = PE(:,M*D1.ca.Order+1:M*(D1.ca.Order+1));
    X.C=C;
    X.D = sqrtm(X.C);
    X.datatype = 'MVAR';
    [S,h,pdc,COH,dtf,DC,pCOH,ddtf,ffdtf, pCOH2, PDCF, coh,GGC,Af,gpdc,GGC2,GDTF,DffDTF,gddtf,ffpdc,gffpdc]=ImaGIN_mvfreqz(X.B,X.A,X.C,D1.ca.Freq,D1.fsample,ErrorCov);
    D1.ca.S=S;
    D1.ca.DTF=dtf;
    D1.ca.ffDTF=ffdtf;
    D1.ca.dDTF=ddtf;
    D1.ca.GdDTF=gddtf;
    D1.ca.PDC=pdc;
    D1.ca.ffPDC=ffpdc;
    D1.ca.GffPDC=gffpdc;
    D1.ca.GPDC=gpdc;
    D1.time=[];
    D1=rmfield(D1,'TimeZero');
else
    D1.ca.S = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.DTF = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.ffDTF = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.dDTF = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.GdDTF = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.PDC = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.ffPDC = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.GffPDC = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.GPDC = zeros(length(D1.ca.channels),length(D1.ca.channels),length(TimeWindow),length(D1.ca.Freq));
    D1.ca.Order=zeros(1,length(TimeWindow));
    for i0=1:length(TimeWindow)
        win=find(time>=TimeWindow(i0)-TimeWindowWidth/2&time<=TimeWindow(i0)+TimeWindowWidth/2);
        
%         %Amplitude normalisation
%         if D1.ca.FlagNormalise
            datatmp=ImaGIN_Normalisation(data(win,:),1);
%         else
%             datatmp=data(win,:);
%         end

        %Model order
%         [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/5));
        [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/8));
%         try
%             [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/4),'fpe');
%         catch
%             [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/6),'fpe');
%         end            
        D1.ca.Order(i0)=size(A,2)/size(A,1);
%         [AR2,RC,PE] = mvar(datatmp,D1.ca.Order(i0));
        AR2=A;
        M = size(AR2,1);
        X.A = [eye(M),-AR2];
        X.B = eye(M);
%         X.C = PE(:,M*D1.ca.Order(i0)+1:M*(D1.ca.Order(i0)+1));
        X.C = C;
        X.D = sqrtm(X.C);
        X.datatype = 'MVAR';
        [S,h,pdc,COH,dtf,DC,pCOH,ddtf,ffdtf, pCOH2, PDCF, coh,GGC,Af,gpdc,GGC2,GDTF,DffDTF,gddtf,ffpdc,gffpdc]=ImaGIN_mvfreqz(X.B,X.A,X.C,D1.ca.Freq,D1.fsample,ErrorCov);
        D1.ca.S(:,:,i0,:)=S;
        D1.ca.DTF(:,:,i0,:)=dtf;
        D1.ca.ffDTF(:,:,i0,:)=ffdtf;
        D1.ca.dDTF(:,:,i0,:)=ddtf;
        D1.ca.GdDTF(:,:,i0,:)=gddtf;
        D1.ca.PDC(:,:,i0,:)=pdc;
        D1.ca.ffPDC(:,:,i0,:)=ffpdc;
        D1.ca.GffPDC(:,:,i0,:)=gffpdc;
        D1.ca.GPDC(:,:,i0,:)=gpdc;
    end
    [tmp t]=min(TimeWindow);
    D1=timeonset(D1,tmp);
end


D1.ca.p=Inf;
n=0;
D1.ca.S_S = [];
D1.ca.DTF_S = [];
D1.ca.ffDTF_S = [];
D1.ca.dDTF_S = [];
D1.ca.GdDTF_S = [];
D1.ca.PDC_S = [];
D1.ca.ffPDC_S = [];
D1.ca.GffPDC_S = [];
D1.ca.GPDC_S = [];
NTimeWindow=max([1 length(TimeWindow)]);
while min(D1.ca.p)>D1.ca.pmax
    n=n+1;
    D1.ca.p=[1:n*NTimeWindow]/(1+n*NTimeWindow);
    %Shift or FT1 surrogates
    switch D1.ca.Surro
        case 'Shift'
            tmp=round(1.9*D1.nsample*(rand(1,D1.nchannels)-0.5));
            Surro=data;
            for i0=1:length(D1.ca.channels)
                Surro(:,i0)=ImaGIN_shift(data(:,i0),tmp(i0),1);
            end
        case 'FT1'
            Surro=data;
            for i0=1:length(D1.ca.channels)
                Surro(:,i0)=ImaGIN_FT1Surrogate(data(:,i0),1);
            end
    end

    if isempty(Baseline)
        Surro=ImaGIN_Normalisation(Surro,1);
        %Find error covariance during the baseline (used for GPDC)
        [w,A,C,SBC,FPE,th]=arfit(Surro,2,ceil(D1.fsample/8));
%         try
%             [w,A,C,SBC,FPE,th]=arfit(Surro,2,ceil(D1.fsample/4),'fpe');
%         catch
%             [w,A,C,SBC,FPE,th]=arfit(Surro,2,ceil(D1.fsample/6),'fpe');
%         end            
%         [AR2,RC,PE] = mvar(Surro,size(A,2)/size(A,1));
%         ErrorCov = PE(:,size(AR2,1)*size(A,2)/size(A,1)+1:size(AR2,1)*(size(A,2)/size(A,1)+1));
        ErrorCov=C;
    else
        Surro=ImaGIN_Normalisation(Surro,1,Baseline);
        %Find error covariance during the baseline (used for GPDC)
        [w,A,C,SBC,FPE,th]=arfit(Surro(Baseline,:),2,ceil(D1.fsample/8));
%         [w,A,C,SBC,FPE,th]=arfit(Surro(Baseline,:),2,ceil(D1.fsample/4),'fpe');
%         try
%             [w,A,C,SBC,FPE,th]=arfit(Surro(Baseline,:),2,ceil(D1.fsample/4),'fpe');
%         catch
%             [w,A,C,SBC,FPE,th]=arfit(Surro(Baseline,:),2,ceil(D1.fsample/6),'fpe');
%         end            
%         [AR2,RC,PE] = mvar(Surro(Baseline,:),size(A,2)/size(A,1));
%         ErrorCov = PE(:,size(AR2,1)*size(A,2)/size(A,1)+1:size(AR2,1)*(size(A,2)/size(A,1)+1));
        ErrorCov = C;
    end

    S_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    DTF_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    ffDTF_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    dDTF_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    GdDTF_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    PDC_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    ffPDC_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    GffPDC_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    GPDC_S = zeros(length(D1.ca.channels),length(D1.ca.channels),NTimeWindow,length(D1.ca.Freq));
    for i0=1:NTimeWindow
        if ~isempty(TimeWindow)
            win=find(time>=TimeWindow(i0)-TimeWindowWidth/2&time<=TimeWindow(i0)+TimeWindowWidth/2);
        else
            win=1:length(time);
        end
        
%         %Amplitude normalisation
%         if D1.ca.FlagNormalise
%             datatmp=ImaGIN_Normalisation(Surro(win,:),1);
%         else
            datatmp=Surro(win,:);
%         end
        
        %Model order
%         [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/5));
        [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/8));
%         try
%             [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/4),'fpe');
%         catch
%             [w,A,C,SBC,FPE,th]=arfit(datatmp,2,ceil(D1.fsample/6),'fpe');
%         end
        Order=size(A,2)/size(A,1);
%         [AR2,RC,PE] = mvar(datatmp,Order);
        AR2=A;
        M = size(AR2,1);
        X.A = [eye(M),-AR2];
        X.B = eye(M);
%         X.C = PE(:,M*Order+1:M*(Order+1));
        X.C=C;
        X.D = sqrtm(X.C);
        X.datatype = 'MVAR';
        [S,h,pdc,COH,dtf,DC,pCOH,ddtf,ffdtf, pCOH2, PDCF, coh,GGC,Af,gpdc,GGC2,GDTF,DffDTF,gddtf,ffpdc,gffpdc]=ImaGIN_mvfreqz(X.B,X.A,X.C,D1.ca.Freq,D1.fsample,ErrorCov);
        S_S(:,:,i0,:)=S;
        DTF_S(:,:,i0,:)=dtf;
        ffDTF_S(:,:,i0,:)=ffdtf;
        dDTF_S(:,:,i0,:)=ddtf;
        GdDTF_S(:,:,i0,:)=gddtf;
        PDC_S(:,:,i0,:)=pdc;
        ffPDC_S(:,:,i0,:)=ffpdc;
        GffPDC_S(:,:,i0,:)=gffpdc;
        GPDC_S(:,:,i0,:)=gpdc;
    end
    D1.ca.S_S=cat(3,D1.ca.S_S,S_S);
    D1.ca.DTF_S=cat(3,D1.ca.DTF_S,DTF_S);
    D1.ca.ffDTF_S=cat(3,D1.ca.ffDTF_S,ffDTF_S);
    D1.ca.dDTF_S=cat(3,D1.ca.dDTF_S,dDTF_S);
    D1.ca.GdDTF_S=cat(3,D1.ca.GdDTF_S,GdDTF_S);
    D1.ca.PDC_S=cat(3,D1.ca.PDC_S,PDC_S);
    D1.ca.ffPDC_S=cat(3,D1.ca.ffPDC_S,ffPDC_S);
    D1.ca.GffPDC_S=cat(3,D1.ca.GffPDC_S,GffPDC_S);
    D1.ca.GPDC_S=cat(3,D1.ca.GPDC_S,GPDC_S);
end

% D1=fsample(D1,1/(TimeWindow(2)-TimeWindow(1)));
% 
D1s=clone(D1,['ar' Pre '_ca1_' fnamedat], [length(D1.ca.channels) size(data,1) D1.ntrials]);
for i1=1:length(D1.ca.chanlabels)
    D1s=chanlabels(D1s,i1,D1.ca.chanlabels{i1});
end
D1s(:,:,k)=data';
clear D2


%  %Récupération des seuils pour chaque couples de channels (16 seuils dont douze significatifs car on ne prends pas en compte la diagonale)
%  %On calcule les seuils pour 3 valeurs différentes de p (4ème dimension de
%  %la matrice résultat.
%  Seuil_DTF=zeros(Nchannels,Nchannels,length(TimeWindow),3);
%  Seuil_dDTF=zeros(Nchannels,Nchannels,length(TimeWindow),3);
%  %Seuils pour la différence
%  Seuil_diff_DTF=zeros(Nchannels,Nchannels,length(TimeWindow),3);
%  Seuil_diff_dDTF=zeros(Nchannels,Nchannels,length(TimeWindow),3);
%  
%  tableau=[0 1 4];    
%  for k=1:length(TimeWindow)
%      a=squeeze(Surrogates_dDTF(:,:,k,:));
%      b=squeeze(Surrogates_DTF(:,:,k,:));
%      c=a-permute(a,[2 1 3]);
%      d=b-permute(b,[2 1 3]);
%      a_sort=sort(real(a),3);
%      b_sort=sort(real(b),3);
%      c_sort=sort(real(c),3);
%      d_sort=sort(real(d),3);
%      for l=1:3
%          a_value=a_sort(:,:,NRealisations-tableau(l));
%          b_value=b_sort(:,:,NRealisations-tableau(l));
%          c_value=c_sort(:,:,NRealisations-tableau(l));
%          d_value=d_sort(:,:,NRealisations-tableau(l));
%          Seuil_dDTF(:,:,k,l)=a_value;
%          Seuil_DTF(:,:,k,l)=b_value;
%          Seuil_diff_dDTF(:,:,k,l)=c_value;
%          Seuil_diff_DTF(:,:,k,l)=d_value;
%      end
%  end
%  
%  
%  D1.ca.Seuil_DTF=Seuil_DTF;
%  D1.ca.Seuil_dDTF=Seuil_dDTF;
%  D1.ca.Seuil_diff_DTF=Seuil_diff_DTF;
%  D1.ca.Seuil_diff_dDTF=Seuil_diff_dDTF;
 
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CM=ConnectivityMatrix(N)

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AC=Crosscorrelation(x,y);

x=x-mean(x);
y=y-mean(y);
N=6;

S=sqrt(sum(x.^2).*sum(y.^2));
NTime=length(x);
AC = zeros(1,2*round(NTime/N)+1);
for i3 =-round(NTime/N):round(NTime/N)
    tmp2=ImaGIN_shift(y,i3);
    AC(i3+round(NTime/N)+1)=sum(x.*tmp2)./S;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==========================================================================
function event = select_events(event, timeseg)
% Utility function to select events according to time segment

if ~isempty(event)
    [time ind] = sort([event(:).time]);

    selectind = ind(time >= timeseg(1) & time <= timeseg(2));

    event = event(selectind);
    for i=1:length(event)
        event(i).time=event(i).time-timeseg(1);
    end
end