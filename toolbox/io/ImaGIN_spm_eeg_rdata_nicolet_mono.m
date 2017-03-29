function D = ImaGIN_spm_eeg_rdata_nicolet_mono(S)
% Converts EEG data from Nicolet .e format to SPM-format with bipolar montage for SEEG
%
% USAGE:   D = ImaGIN_spm_eeg_rdata_nicolet_mono(S)

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


D = [];


try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.e$', 'Select Nicolet .e file');
end
cd(spm_str_manip(Fdata,'h'))
FdataName=spm_str_manip(Fdata,'t');
FdataName=FdataName(1:end-2);

try
    FileOut=S.FileOut;
catch
    FileOut=spm_input('Name of converted file', '+1', 's');
end

% which species?
try
    S.Atlas;
catch
    Ctype = {
            'Human',...
            'Rat',...
            'Mouse'};
    str   = 'Select atlas';
	Sel   = spm_input(str, '+1', 'm', Ctype);
	S.Atlas = Ctype{Sel};
end


%Read Nicolet file
obj=NicoletFile(Fdata);

%find channels with EEG and annotations
try
    S.channel;
catch
    S.channel = str2num(spm_input('Index of channels to read', '+1', 's'));
end
if isempty(S.channel)  %for the moment reads only electrodes with -, assuming it's bipolar montage
    Annotations=[];
    SEEGBip=[];
    for i1=1:length(obj.segments.chName)
        if strcmp(obj.segments.chName{i1},'Annotations')
            Annotations=i1;
        end
        if ~isempty(strfind(obj.segments.chName{i1},'-'))
            SEEGBip=[SEEGBip i1];
        end
    end
    S.channel=SEEGBip;
end

%Save name of electrodes
for i1=1:length(S.channel)
    D.channels.label{i1,1}=obj.segments.chName{S.channel(i1)};
end
          
D.channels.label=lower(D.channels.label);

S.reference=[];

try
    S.coarse;
catch
%     S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
    spm_input(sprintf('Sampling rate of acquisition = %d Hz',obj.segments.samplingRate(S.channel(1))),'+1','d');
    tmp=str2num(spm_input('Sampling rate of conversion [Hz]', '+1', 's',obj.segments.samplingRate(S.channel(1))));
    S.coarse=max([ceil(obj.segments.samplingRate(S.channel(1))/tmp) 1]);
    spm_input(sprintf('Effective sampling rate of conversion = %d Hz',obj.segments.samplingRate(S.channel(1))/S.coarse),'+1','d');
end

spm('Pointer','Watch'); drawnow;

D.Atlas=S.Atlas;


F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');

% Read data
if ~isempty(S.channel)
    Data=zeros(length(S.channel),obj.segments.samplingRate(S.channel(1))*obj.segments.duration);
    for i1=1:length(S.channel)
        if obj.segments.samplingRate(S.channel(i1)) == obj.segments.samplingRate(S.channel(1))
           Data(i1,:) = getdata(obj, 1, [1 obj.segments.samplingRate(S.channel(i1))*obj.segments.duration], S.channel(i1))';
        else 
            d=getdata(obj, 1, [1 obj.segments.samplingRate(S.channel(i1))*obj.segments.duration], S.channel(i1))';
            Data(i1,:) = cat(2,d,zeros(1,length(Data(1,:))-length(d)));
        end
    end
else
    Data=[];
    ok=1;
    i1=0;
    S.channel=[];
    while ok
        i1=i1+1;
        try
            Data(i1,:) = getdata(obj, 1, [1 obj.segments.samplingRate(i1)*obj.segments.duration], i1)';
            S.channel=[S.channel i1];
        catch
            ok=0;
        end
    end
end



% Read data and time strings
D.descrip.date = obj.segments.startDate;
D.descrip.time = obj.segments.startTime;

% Read number of channels
Nchannels = size(Data,1);

D.timeOnset=0;

% Read AD rate
D.Fsample=obj.segments.samplingRate(S.channel(1));

%DownSampling
if S.coarse>1
    for i1=1:size(Data,1)
        Data(i1,:)= ImaGIN_lowpassFilter(Data(i1,:),D.Fsample,D.Fsample/(2*S.coarse)-2);    %antialiasing,
    end
end
Data=Data(:,1:S.coarse:end);
D.Fsample=D.Fsample/S.coarse;

time=[0:size(Data,2)-1]./D.Fsample;
D.time=time;


% Number of expected samples
Nsamples = size(Data,2); 

% store only relative name of channel template file
D.channels = repmat(struct('bad', 0,'label',D.channels.label,'type','EEG'), 1, 1);

%rename channels for Bucharest
for i1=1:length(D.channels)
    a_string=D.channels(i1).label;
    tmp=find(a_string=='-');
    if ~isempty(tmp)
        tmp1=a_string(1:tmp-1);
        tmp2=a_string(tmp+1:end);
        tmp=find(tmp1=='0');
        if ~isempty(tmp)
            if tmp<length(tmp1)
                tmp1=tmp1(setdiff(1:length(tmp1),tmp));
            end
        end
        tmp=find(tmp2=='0');
        if ~isempty(tmp)
            if tmp<length(tmp2)
                tmp2=tmp2(setdiff(1:length(tmp2),tmp));
            end
        end
        D.channels(i1).label=[tmp2 tmp1];
    end
end


D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
nsampl=Nsamples;
nchan=D.Nchannels;

%Events
D.trials.onset = 1/D.Fsample;


%Continuous
D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
% physically initialise file
D.data(end,end) = 0;
offset = 1;
nblocksamples = size(Data,2);
D.data(:, offset:(offset+nblocksamples-1)) = full(Data);


% %Electrodes

for i1=1:length(D.channels)
    D.sensors.eeg.label{1,i1}=D.channels(i1).label;
end

D.sensors.eeg.unit='mm';
D.sensors.eeg.pnt=NaN*zeros(length(D.channels),3);
D.sensors.eeg.type='eeg';
D.sensors.chantype='eeg';

D.Atlas=S.Atlas;

%--------- Create meeg object
[DirOut,NameOut,~]=fileparts(FileOut);
D.fname = [NameOut '.mat'];
D.path = DirOut;

D = meeg(D);
save(D)

spm('Pointer','Arrow');






