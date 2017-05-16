function ImaGIN_spike_compute(S)
% Compute mean firing rate
%
% USAGE:   D = ImaGIN_spike_compute(S)

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


%% ===== LOAD INPUT FILE =====
% Get input file
try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select EEG mat file');
end
% Change mouse pointer
spm('Pointer', 'Watch'); drawnow;
% Load meeg object
Dmat = load(DD);
% If the spike field is not available, try to load it with the "load" function
if ~isfield(Dmat, 'D') || ~isfield(Dmat.D, 'spike')
    error('No spike data');
end
% Copy spike field
Dspike = Dmat.D.spike;
% Convert to MEEG object
D = meeg(Dmat.D);


%% ===== OPTIONS =====
%Recalculate seizure
NElectrode=length(Dspike.name);
try
    SelChan = S.Channel;
catch
    SelChan = spm_input('Select channel(s) ', 1, 'i',1:NElectrode);
end

try
    SelNeuro = S.Neurone;
catch
    for i1=1:length(SelChan)
        NNeurone=max(Dspike.markers{SelChan(i1)}(:,1,1,1));
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



%% ===== COMPUTATION =====
%Save parameter analysis
for i1=1:length(SelChan)
    Dspike.Analysis{SelChan(i1)}.SelNeuro=SelNeuro;
    Dspike.Analysis{SelChan(i1)}.TimeWindow=TimeWindow;
    Dspike.Analysis{SelChan(i1)}.TimeWindowWidth=TimeWindowWidth;
    if isempty(TimeWindow)
        Dspike.Analysis{SelChan(i1)}.TimeResolution=TimeResolution;
        Dspike.Analysis{SelChan(i1)}.Start=Start;
        Dspike.Analysis{SelChan(i1)}.End=End;
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
    Index = find(Dspike.timings{SelChan(i1)}>=Start & Dspike.timings{SelChan(i1)}<=End);
    timings = Dspike.timings{SelChan(i1)}(Index);
    markers = double(Dspike.markers{SelChan(i1)}(Index,1,1,1));



    if Instantaneous
        %Instantaneous firing rate for each neurone
        for i3=1:length(SelNeuro{i1})
            n=n+1;
            D.Firing.Name{n}=[Dspike.name{i1} ' - ' num2str(SelNeuro{i1}(i3))];
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
            D.Firing.Name{n}=[Dspike.name{i1} ' - ' num2str(SelNeuro{i1}(i3))];
        end
    end
end
if Instantaneous
    S=zeros(size(Data,1),length(Time));
    for i1=1:size(Data,1)
        S(i1,:) = interp1(min(Time):TimeBinDuration:max(Time),Data(i1,:),Time,'linear');
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


%% ===== SAVE FILE =====
% Get file name
Dfile = fullfile(D);
% Copy results of spike analysis in the file structure
Dmat.D.Firing.data = S;
D = Dmat;
% Save file
save(Dfile, 'D');

% Restore mouse pointer
spm('Pointer', 'Arrow');



