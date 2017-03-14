function D = ImaGIN_spm_eeg_rdata_ascii(S)
% converts EEG data from ascii - to SPM-format
% FORMAT D = spm_eeg_rdata_ascii(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of Micromed-file
% Fchannels   - filename of channel template
% reference   - name of reference channel(s)
% channel     - index of channel(s) to read
%_______________________________________________________________________
% 
% ImaGIN_spm_eeg_rdata_micromedsys3 reads a continuous *.trc Micromed file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rdata.m 317 2005-11-28 18:31:24Z stefan $


try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.asc$', 'Select ASCII file');
end
cd(spm_str_manip(Fdata,'h'))

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

try 
    CreateTemplate=S.CreateTemplate;
catch
    CreateTemplate = spm_input('Create template ','+1','Yes|No');
end

switch CreateTemplate
    case{'No'}
        try
            Fchannels = S.Fchannels;
        catch
            Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
        end
    case{'Yes'}
        P=spm_str_manip(Fdata,'h');
        try
            filename=S.filenameName;
        catch
            filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
        end
        fid=fopen(filename,'r');
        n=0;
        while 1
            tmp=fgetl(fid);
            if ~ischar(tmp)
                break
            else
                n=n+1;
                Cnames{n,1}=tmp;
            end
        end
        fclose(fid);
        try
            filename=S.filenamePos;
        catch
            filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
        end
        Position=load(filename);

        Cpos2=Position';
        Cpos=Position(:,1:2)';
        Cpos=Cpos-min(Cpos(:));
        Cpos=Cpos./max(Cpos(:));

        try
            Species=S.Atlas;
        catch
            Species = spm_input('Species ','+1','Human|Rat|Mouse');
        end
        switch Species
            case{'Rat','Mouse'}
                Cpos2=10*Cpos2;
        end

        Files=what(fullfile(spm('dir'),'EEGTemplates'));
        ok=1;
        while ok
            try
                Filename=S.MontageName;
                ok=0;
            catch
                Filename=spm_input('Name of montage', '+1','s');
                ok=0;
                for i1=1:length(Files.mat)
                    if strcmp(lower(Files.mat{i1}),lower([Filename '.mat']))
                        spm_input('Specified file already exists', '+1','d')
                        ok=1;
                        break
                    end
                end
            end
        end
        Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
        Rxy=1;
        Nchannels=length(Cnames);
        if str2num(version('-release'))>=14
            save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')
        else
            save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species')
        end

        Fchannels=Filename;
        Csetup = load(Fchannels);
end


% % There doesn't seem to be information in the CNT-file which channel(s) is the reference
% % so ask for it now
% try
%     S.reference;
% catch
%     S.reference = spm_input('Input reference channel name', '+1', 's');
% end
S.reference=[];

try
    Bipolar=S.Bipolar;
catch
    Bipolar = spm_input('Bipolar derivations? ','+1','Yes|No');
end
switch Bipolar
    case{'No'}
        try
            S.channel;
        catch
            S.channel = str2num(spm_input('Index of channels to read', '+1', 's'));
        end
    case{'Yes'}
        try
            S.bipole;
        catch
            P=spm_str_manip(Fdata,'h');
            tmp = spm_select(1, '\.txt$', 'Select txt file for bipolar derivations', {}, P);
            S.bipole=load(tmp);
        end
end

try
    S.coarse;
catch
    S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
end

spm('Pointer','Watch'); drawnow;

D = [];
D.Atlas=S.Atlas;

if strcmp(CreateTemplate,'No')
    Csetup = load(Fchannels);
%     if length(S.channel)~=Csetup.Nchannels&~isempty(S.channel)
%         error('Number of read channels does not match the template')
%     end
end


F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');


% file size
Dn = dir(Fdata);
Sf = Dn.bytes;

% Read data
Data=load(Fdata);
% if size(Data,2)<size(Data,1)
%     Data=Data';
% end
% if ~strcmp(CreateTemplate,'Yes')
%     S.channel=1:size(Data,1);
% end

if isfield(S,'channel')
    Data=Data(S.channel,:);
elseif isfield(S,'bipole')
    Data=Data(S.bipole(:,1),:)-Data(S.bipole(:,2),:);
end

% if strcmp(CreateTemplate,'Yes')
%     P=spm_str_manip(Fdata,'h');
%     filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
%     fid=fopen(filename,'r');
%     n=0;
%     while 1
%         tmp=fgetl(fid);
%         if ~ischar(tmp)
%             break
%         else
%             n=n+1;
%             Cnames{n,1}=tmp;
%         end
%     end
%     fclose(fid);
%     filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
%     Position=load(filename);
%     
%     Cpos2=Position';
%     Cpos=Position(:,1:2)';
%     Cpos=Cpos+min(Cpos(:));
%     Cpos=Cpos./max(Cpos(:));
%     
%     Species = spm_input('Create template ','+1','Rat|Mouse|Human');
%     
%     Files=what(fullfile(spm('dir'),'EEGTemplates'));
%     ok=1;
%     while ok
%         Filename=spm_input('Name of montage', '+1','s');
%         ok=0;
%         for i1=1:length(Files.mat)
%             if strcmp(lower(Files.mat{i1}),lower([Filename '.mat']))
%                 spm_input('Specified file already exists', '+1','d')
%                 ok=1;
%                 break
%             end
%         end
%     end
%     Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
%     Rxy=1;
%     save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')
% 
%     Fchannels=Filename;
%     Csetup = load(Fchannels); 
% end

% Read data and time strings
D.descrip.date = ' / / ';
D.descrip.time = ' / / ';

% Read number of channels
Nchannels = size(Data,1);

%Time zero (OD)
D.TimeZero=1;

% Read AD rate
try
    D.Fsample=S.Radc;
catch
    D.Fsample = str2num(spm_input('Sampling frequency [Hz] ', '+1', 's'));
end


%DownSampling
if S.coarse>1
    for i1=1:size(Data,1)
        Data(i1,:)= ImaGIN_lowpassFilter(Data(i1,:),D.Fsample,D.Fsample/(2*S.coarse)-2);    %antialiasing,
    end
end
Data=Data(:,1:S.coarse:end);
D.Fsample=D.Fsample/S.coarse;
time=[0:size(Data,2)-1]./D.Fsample;

%Events
try 
    Nevent=S.Nevent;
catch
    Nevent = str2num(spm_input('Number of event types', '+1', 's',0));
end
Nevents=0;
Ec=[];
Et=[];
for i1=1:Nevent
    try
        D.trials.label{i1}=S.event(i1); %I changed S.event{i1}
    catch
        D.trials.label{i1}=spm_input(sprintf('Name of event %d',i1), '+1', 's');
    end
        
    try
        tmp=deblank(S.event_file(i1,:));
    catch
        tmp = spm_select(1, '\.txt$', sprintf('Select txt file with the timing (sec) of event %d',i1));
    end
    
    tmp=load(tmp);
    Nevents=Nevents+length(tmp);
    Ec=[Ec;i1*ones(length(tmp),1)];
    for i2=1:length(tmp)
        Et=[Et;unique(find(abs(time-tmp(i2))==min(abs(time-tmp(i2)))))];
    end
end
    
% Number of expected samples
Nsamples = size(Data,2); 

% store only relative name of channel template file
D.channels = repmat(struct('bad', 0), 1, 1);


% Read the electrode names
if isfield(S,'channel')
    for i2=1:length(S.channel)
        D.channels(i2).label=Csetup.Cnames(i2);
    end
else
    for i2=1:length(Csetup.Cnames)
        D.channels(i2).label=Csetup.Cnames(i2);
    end
end


% Map name of channels to channel order specified in channel template file
for i = 1:length(D.channels)
	index = [];
	for j = 1:Csetup.Nchannels
		if ~isempty(find(strcmpi(deblank(D.channels(i).label), Csetup.Cnames{j})))
			index = [index j];
		end
	end
    
	if isempty(index)
		warning(sprintf('No channel named %s found in channel template file.', D.channels(i).label));
	else
		% take only the first found channel descriptor
		D.channelOrder(i) = index(1);
    end
end

% D.data.fnamedat = [S.SaveFile '.dat'];
% D.data.datatype = 'float32-le';

D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
nsampl=Nsamples;
nchan=D.Nchannels;


if ~isempty(Ec)
    for i1=1:max(Ec)
        D.trials(i1).events=find(Ec==i1);
        D.trials(i1).onset=Etsec(find(Ec==i1));
    end
else
    D.trials.label = 'Undefined';
    D.trials.events = [];
    D.trials.onset = 1/D.Fsample;
end

%Continuous
D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
% physically initialise file
D.data(end,end) = 0;
offset = 1;
nblocksamples = size(Data,2);
D.data(:, offset:(offset+nblocksamples-1)) = full(Data);


%Electrodes
try
    D.sensors.eeg.pnt=Csetup.Cpos2(:,ChannelOrder)';
catch
    D.sensors.eeg.pnt=Csetup.Cpos2(:,:)';
end
for i1=1:length(D.channels)
%     D.sensors.eeg.label{1,i1}=D.channels(i1).label;
    D.sensors.eeg.label(i1)=D.channels(i1).label;
end
D.sensors.eeg.unit='mm';
% D.sensors.eeg.type='eeg';

D.Atlas=S.Atlas;

%--------- Create meeg object
[DirOut,NameOut,~]=fileparts(FileOut);
D.fname = [NameOut '.mat'];
D.path = DirOut;

D = meeg(D);
save(D)

spm('Pointer','Arrow');
