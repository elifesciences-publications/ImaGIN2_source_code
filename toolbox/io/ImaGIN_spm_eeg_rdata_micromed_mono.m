function D = ImaGIN_spm_eeg_rdata_micromed_mono(S)
% converts EEG data from Micromed - to SPM-format
% FORMAT D = spm_eeg_rdata_micromed(S)
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
	Fdata = spm_select(1, '\.TRC$', 'Select Micromed file');
end
cd(spm_str_manip(Fdata,'h'))

% try
%     System = S.System;
% catch
%     Ctype = {
%             'System 3',...
%             'System 4',...
%             'System 98'};
%     str   = 'Select Micromed system';
% 	Sel   = spm_input(str,'+1', 'm', Ctype);
% 	System = Ctype{Sel};
% %     System = spm_input('Select Micromed system','+1','System 3|System 4|System 98'),
% end

% which species?
try
    S.Atlas;
catch
%     Ctype = {
%             'Human',...
%             'Rat',...
%             'Mouse'};
    str   = 'Select atlas';
% 	Sel   = spm_input(str, '+1', 'm', Ctype);
% 	S.Atlas = Ctype{Sel};
    S.Atlas=spm_input(str, '+1','Human|Rat|Mouse');
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
    S.channel;
catch
    S.channel = str2num(spm_input('Index of channels to read', '+1', 's'));
end

try
    S.coarse;
catch
    S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
end

try
    FileOut=S.FileOut;
catch
    FileOut=spm_input('Name of converted file', '+1', 's');
end


spm('Pointer','Watch'); drawnow;


%--------- Start making the header
D = [];

F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');


% fp = fopen(Fdata, 'r');

% file size
Dn = dir(Fdata);
Sf = Dn.bytes;

% % Read data
% switch System
% %     case{'System 3','System 4'}
    try
%         [Data,fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=readtrc(Fdata,System,S);
        [Data,fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=readtrc(Fdata,S);
        try
            tmp=Data(v_neworder_rev,:);
            Data=tmp;
            if S.coarse>1
                for i1=1:size(Data,1)
                    %                 Data(i1,:)= ImaGIN_lowpassFilter(Data(i1,:),fff_raw.file_format.sampling_rate,fff_raw.file_format.sampling_rate/(2*S.coarse)-2);    %antialiasing,
                    Data(i1,:)= ImaGIN_bandpassFilter(Data(i1,:),fff_raw.file_format.sampling_rate,0.6,fff_raw.file_format.sampling_rate/(2*S.coarse)-2);    %antialiasing,
                end
                Data=Data(:,1:S.coarse:end);
            end
        catch
            if S.coarse~=1
                tmp=Data(:,1:S.coarse:end);
                Data=tmp;
            end
            tmp=Data(v_neworder_rev,:);
            Data=tmp;
            clear tmp
        end            
        time=inv(fff_raw.file_format.sampling_rate)*[0:size(Data,2)-1];
        if isempty(S.channel)
            S.channel=1:size(Data,1);
        end
        Data=Data(S.channel,:);

    
        try
            loadevents=S.loadevents;
        catch
            loadevents=spm_input('Load notes from .TRC? ','+1','yes|no');
        end
        if strcmp(loadevents,'yes')
            tmploc=[];
            for i1=1:length(ss_note)
                if ss_note(i1).time_in_sample>0
                    if ~(strcmp(ss_note(i1).text(1),'*'))
                        tmploc=[tmploc i1];
                    end
                end
            end
            ss_note=ss_note(tmploc);
%             D.event.ss_note=ss_note;
            Event=zeros(1,length(ss_note));
            Et=zeros(1,length(ss_note));
            for i1=1:length(Event)
                %             Event(i1)=str2num(TRC.event(i1).type);
                Event(i1)=i1;
                Et(i1)=max([ss_note(i1).time_in_sample 1]);
                D.events.name{i1}=deblank(ss_note(i1).text);
            end
            Event1=unique(Event);
            Ec=Event;
            %         for i1=1:length(Event1)
            %             tmp=find(Event==Event1(i1));
            %             Ec(tmp)=i1;
            %             if length(Event1)<10
            %                 try
            %                     D.events.name{i1}=S.event{i1};
            %                 catch
            %                     D.events.name{i1}=spm_input(sprintf('Name of event %d',i1), '+1', 's');
            %                 end
            %             else
            %                 D.events.name{i1}='';
            %             end
            %         end
        else
            Et=[];
            Ec=[];
            
        end
    
%     case{'System 98'}
    catch
        clear Param
        Param.filename=Fdata;
        Param.chans=num2str(S.channel);
        Param.chan_adjust_status=0;
        Param.chan_adjust='';
        try
            Param.loadevents.state=S.loadevents;
        catch
            Param.loadevents.state=spm_input('Load events from .TRC? ','+1','yes|no');
        end
        switch Param.loadevents.state
            case{'yes'}
                try
                    tmp=S.eventtype;
                catch
                    tmp=spm_input('Type of trigger? ','+1','Analog|Digital|Both');
                end
                switch tmp
                    case{'Analog'}
                        Param.loadevents.type='marker';
                        Param.loadevents.dig_ch1='';
                        Param.loadevents.dig_ch2='';
                        Param.loadevents.dig_ch1_label='';
                        Param.loadevents.dig_ch2_label='';
                    case{'Digital'}
                        Param.loadevents.type='eegchan';
                        Param.loadevents.dig_ch1=spm_input('Channel name for marker 1 ', '+1', 's');
                        Param.loadevents.dig_ch1_label=spm_input('Name for marker 1 ', '+1', 's',1);
                        Param.loadevents.dig_ch2=spm_input('Channel name for marker 2 ', '+1', 's');
                        Param.loadevents.dig_ch2_label=spm_input('Name for marker 2 ', '+1', 's',2);
                    case{'Both'}
                        Param.loadevents.type='both';
                end                        
            case{'no'}
                Param.loadevents.type='';
                Param.loadevents.dig_ch1='';
                Param.loadevents.dig_ch2='';
                Param.loadevents.dig_ch1_label='';
                Param.loadevents.dig_ch2_label='';
        end
        TRC=readtrc98(Param);
        Data=TRC.data;
        if isfield(S,'bipole')
            DataTmp=zeros(size(Data,1)/2,size(Data,2));
%             for i1=1:size(DataTmp,1)
%                 DataTmp(i1,:)=Data(i1,:)-Data(i1+size(Data,1)/2,:);
%             end
            n=0;
            for i1=1:2:2*size(DataTmp,1)
                n=n+1;
                if S.channel(i1)~=S.channel(i1+1)
                    DataTmp(n,:)=Data(i1,:)-Data(i1+1,:);
                else
                    DataTmp(n,:)=Data(i1,:);
                end
            end
            Data=DataTmp;
        end
        if S.coarse>1
            for i1=1:size(Data,1)
                Data(i1,:)= ImaGIN_lowpassFilter(Data(i1,:),TRC.srate,TRC.srate/(2*S.coarse)-2);    %antialiasing,
            end
            Data=Data(:,1:S.coarse:end);
        end
        Event=zeros(1,length(TRC.event));
        Et=zeros(1,length(TRC.event));
        for i1=1:length(Event)
            Event(i1)=str2num(TRC.event(i1).type);
            Et(i1)=TRC.event(i1).latency;
        end
        Event1=unique(Event);
        Ec=Event;
        for i1=1:length(Event1)
            tmp=find(Event==Event1(i1));
            Ec(tmp)=i1;
            if length(Event1)<10
                try
                    D.events.name{i1}=S.event{i1};
                catch
                    D.events.name{i1}=spm_input(sprintf('Name of event %d',i1), '+1', 's');
                end
            else
                D.events.name{i1}='';
            end                
        end
    end
% end

% % Read data and time strings
% switch System
%     case{'System 3','System 4'}
%         D.descrip.date = [num2str(fff_raw.recording_date.year) '/' num2str(fff_raw.recording_date.month) '/' num2str(fff_raw.recording_date.day)];
%         D.descrip.time = ' / / ';
%         D.Fsample = fff_raw.file_format.sampling_rate;
% 
%     case{'System 98'}
%         D.descrip.date = ' / / ';
%         D.descrip.time = ' / / ';
%         D.Fsample = TRC.srate;
% end
try
    D.descrip.date = [num2str(fff_raw.recording_date.year) '/' num2str(fff_raw.recording_date.month) '/' num2str(fff_raw.recording_date.day)];
    D.descrip.time = ' / / ';
    D.Fsample = fff_raw.file_format.sampling_rate;
catch
    D.descrip.date = ' / / ';
    D.descrip.time = ' / / ';
    D.Fsample = TRC.srate;
end


% Read number of channels
Nchannels=size(Data,1);
nchan=Nchannels;

% Data=Data(:,1:S.coarse:end);
D.Fsample=D.Fsample/S.coarse;

D.timeOnset=0;

% %Events
% switch System
%     case{'System 3','System 4'}
%         
%         if strcmp(loadevents,'no')
%             try
%                 Nevent=S.NeventType;
%             catch
%                 Nevent = str2num(spm_input('Number of event types', '+1', 's',0));
%             end
%             
%             Nevents=0;
%             Ec=[];
%             Et=[];
%             for i1=1:Nevent
%                 try
%                     D.trials.label{i1}=S.event(i1); %I changed S.event{i1}
%                 catch
%                     D.trials.label{i1}=spm_input(sprintf('Name of event %d',i1), '+1', 's');
%                 end
%                 
%                 try
%                     tmp=deblank(S.event_file(i1,:));
%                 catch
%                     tmp = spm_select(1, '\.txt$', sprintf('Select txt file with the timing (sec) of event %d',i1));
%                 end
%                 
%                 tmp=load(tmp);
%                 Nevents=Nevents+length(tmp);
%                 Ec=[Ec;i1*ones(length(tmp),1)];
%                 for i2=1:length(tmp)
%                     Et=[Et;unique(find(abs(time-tmp(i2))==min(abs(time-tmp(i2)))))];
%                 end
%             end
%         end
% end

Et=round(Et/S.coarse);

% Number of expected samples
Nsamples = size(Data,2); 

if isfield(D,'time')
    time=D.time;
else
    time=[0:1/D.Fsample:(Nsamples-1)/D.Fsample];
    time=time-D.timeOnset;
end
Etsec=time(Et);
% D.time=time;

D.channels = repmat(struct('bad', 0), Nchannels,1);


% Read the electrode names
% switch System
%     case{'System 3','System 4'}
%         Name={};
%         for i1=1:length(v_neworder_rev)
%             
%             Name{i1,1}=ss_channel.values(v_neworder_rev(i1)).value;
% 
%             
%             %Lyon
%             Name{i1,1}=lower(Name{i1,1}(~isspace(Name{i1,1})));
%             
%             %Rennes
%             Name{i1,1}=Name{i1,1}(setdiff(1:length(Name{i1,1}),strfind(Name{i1,1},'.')));
%             
%         end
%         for i1=1:length(S.channel)
%             D.channels(i1).label=Name(S.channel(i1));
%             switch D.channels(i1).label{1}
%                 case{'ecg1','ecg2','ekg1','ekg2'}
%                     D.channels(i1).type='ECG';
%                 otherwise
%                     D.channels(i1).type='EEG';
%             end
%         end
%     case{'System 98'}
%         for i2=1:length(Csetup.Cnames)
%             D.channels(i2).label=Csetup.Cnames(i2);
%             D.channels(i2).type='EEG';
%         end
% end
try
    Name={};
    for i1=1:length(v_neworder_rev)
        
        Name{i1,1}=deblank(ss_channel.values(v_neworder_rev(i1)).value);
        
        
        %Lyon
        Name{i1,1}=lower(Name{i1,1}(~isspace(Name{i1,1})));
        
        %Rennes
        Name{i1,1}=Name{i1,1}(setdiff(1:length(Name{i1,1}),strfind(Name{i1,1},'.')));
        
        
        
        %remove zeros
        a=Name{i1,1};
        v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
        if length(v_num)>1
            if v_num(2)==v_num(1)+1 && strcmp(a(v_num(1)),'0')
                Name{i1,1}=a(setdiff(1:length(a),v_num(1)));
            end
        end

    end
    isecg=[];
    for i1=1:length(S.channel)
        D.channels(i1).label=Name(S.channel(i1));
        switch D.channels(i1).label{1}
            case{'ecg1','ecg2','ekg1','ekg2'}
                D.channels(i1).type='ECG';
                isecg(end+1)=i1;
            otherwise
                D.channels(i1).type='EEG';
        end
    end
catch
    for i2=1:length(Csetup.Cnames)
        D.channels(i2).label=Csetup.Cnames(i2);
        D.channels(i2).type='EEG';
    end
end

%reorder to put ECG channels at the end
reorder=[setdiff(1:length(S.channel),isecg) isecg];
D.channels=D.channels(reorder);
Data=Data(reorder,:);

D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
nsampl=Nsamples;


%Assume continuous data
D.trials.label = 'Undefined';
D.trials.events = [];
D.trials.onset = 1/D.Fsample;

if ~isempty(Ec)
    evt = [];
    for i0=1:length(Et)
        evt(i0).type  = D.events.name{Ec(i0)};
        evt(i0).time = time(Et(i0));
        evt(i0).value= Et(i0);
    end
end

    
%Continuous
D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
% physically initialise file
D.data(end,end) = 0;
offset = 1;
nblocksamples = size(Data,2);
D.data(:, offset:(offset+nblocksamples-1)) = full(Data);


%Electrodes
for i1=1:length(D.channels)
%     D.sensors.eeg.label{1,i1}=D.channels(i1).label;
    D.sensors.eeg.label(i1)=D.channels(i1).label;
end
D.sensors.eeg.unit='mm';
D.sensors.eeg.pnt=NaN*zeros(length(D.channels),3);
D.sensors.eeg.type='eeg';
D.sensors.chantype='eeg';

D.Atlas=S.Atlas;

%--------- Create meeg object
[DirOut,NameOut,~]=fileparts(FileOut);
D.fname = [NameOut '.mat'];
D.path=DirOut;
D.type='continuous';

D = meeg(D);
try
    D = events(D, 1, evt);
end
% for i=1:length(D.chanlabels)
%     switch D.chanlabels{i}{1}
%         case{'ecg1','ecg2','ekg1','ekg2'}
%             D=chantype(D,i,'ECG');
%     end
% end


save(D);

spm('Pointer','Arrow');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [data,fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=readtrc(filename,System,S);
function [data,fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=readtrc(filename,S);

% switch System
%     case {'System 3'}
%         [fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=rdhd_sys3(filename);
%     case {'System 4'}
%         [fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=rdhd_sys4(filename);
% end
FileOut=S.FileOut;
try
    [fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=rdhd_sys4(filename, FileOut);
catch
    [fff_raw, ss_channel, s_adress_of_data,v_neworder, v_neworder_rev,m_bipole,ss_note]=rdhd_sys3(filename, FileOut);
end
s_chan=fff_raw.file_format.number_of_channels;
s_chan=1;

f_trc=fopen(filename,'r','n','windows-1252');
sc_an=1;
fseek(f_trc,s_adress_of_data+(s_chan-1)*fff_raw.file_format.size_of_data_bytes,-1)

% data=fread(f_trc,[1,inf],'int16',fff_raw.file_format.size_of_data_bytes*(fff_raw.file_format.number_of_channels-1));

% switch System
%     case {'System 3'}
%         data=fread(f_trc,[fff_raw.file_format.number_of_channels,inf],'int16');
%     case {'System 4'}
%         data=fread(f_trc,[fff_raw.file_format.number_of_channels,inf],'uint16');
% end

data=fread(f_trc,[fff_raw.file_format.number_of_channels,inf],'uint16');

fclose(f_trc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fff_raw, ss_channel, s_adress_of_data, v_neworder, v_neworder_rev, m_bipole, ss_note]=rdhd_sys3(a_trcfilename, FileOut)

% This function reads the header of a micromed .trc file provided that this
% header is of type 4
% It creates a corresponding hdr file, (ascii) with all this info

% exple: [fff_raw, ss_channel, s_adress_of_data,m_bipole]=rdhd_sys4('d:/data2005/eeg_40.trc');

% returns a structure fff_raw that has the characteristics of a fff_raw ella data to describe the exam
% in this directory, we have created a directory a_dirname with the name of the patient and the data of recording
% that's where this hdr file will go.

% get file name and path

[a_path,a_name,a_ext]=fileparts(a_trcfilename);
a_dirname=a_path;

% get the file size
a_com=['ss_d=dir(''' a_trcfilename ''');'];
eval(a_com);
file_size=ss_d.bytes;

% open file
f_trc=fopen(a_trcfilename,'r','n','windows-1252');
% a_hdrname=strrep(lower(a_trcfilename),'.trc','.hdr');
a_trcfilename(end-2:end)=lower(a_trcfilename(end-2:end));
a_hdrname=[FileOut '.hdr'];
f_hdr=fopen(a_hdrname,'w');

fprintf(f_hdr,'%s\n','[File_Size]');
a_string=['File_Size_in_Bytes=' num2str(file_size)];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get title
fprintf(f_hdr,'%s\n','[Title]');
a_s=fread(f_trc,32,'char');
a_string=['Title=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get lab
fprintf(f_hdr,'%s\n','[Lab]');
a_s=fread(f_trc,32,'char');
a_string=['Lab=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get patient info
fprintf(f_hdr,'%s\n','[Patient]');
a_s=fread(f_trc,22,'char');
a_string=['Last_Name=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.name=deblank(char(a_s'));

a_s=fread(f_trc,20,'char');
a_string=['First_Name=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.surname=deblank(char(a_s'));

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Month=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Day=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Year=' num2str(1900+a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.yob=(1900+a_s);

a_s=fread(f_trc,19,'char');
fprintf(f_hdr,'%s\n','');

fff_raw.subject.handedness='?'; % not in the header
fff_raw.subject.sex='?'; % not in the header

% recording date
fprintf(f_hdr,'%s\n','[Recording Date and Time]');

a_s1=fread(f_trc,1,'uchar');
a_string=['Recording_Month=' num2str(a_s1) ];
fprintf(f_hdr,'%s\n',a_string);

a_s2=fread(f_trc,1,'uchar');
a_string=['Recording_Day=' num2str(a_s2) ];
fprintf(f_hdr,'%s\n',a_string);

a_s3=fread(f_trc,1,'uchar');
a_string=['Recording_Year=' num2str(1900+a_s3) ];
fprintf(f_hdr,'%s\n',a_string);

fff_raw.recording_date.year=1900+a_s3;
fff_raw.recording_date.month=a_s1;
fff_raw.recording_date.day=a_s2;

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Hour=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Min=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Sec=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.recording_site.lab='Grenoble';
fff_raw.recording_site.pi='JP_Lachaux';
a_sname=fff_raw.subject.name;
a_sname=a_sname(1:min(3,length(a_sname)));
fff_raw.exam=['exp_' a_sname '_'  int2str(fff_raw.recording_date.day) '_' int2str(fff_raw.recording_date.month) '_' int2str(fff_raw.recording_date.year)];


% get the acquisition equipment
ss_ascii(1).txt='BQ124 - 24 channels headbox, Internal Interface';
ss_ascii(2).txt='MS40 - Holter recorder';
ss_ascii(3).txt='BQ132S - 32 channels headbox, Internal Interface';
ss_ascii(7).txt='BQ124 - 24 channels headbox, BQ CARD Interface';
ss_ascii(8).txt='SAM32 - 32 channels headbox, BQ CARD Interface';
ss_ascii(9).txt='SAM25 - 25 channels headbox, BQ CARD Interface';
ss_ascii(10).txt='BQ132S R - 32 channels reverse headbox, Internal Interface';
ss_ascii(11).txt='SAM32 R - 32 channels reverse headbox, BQ CARD Interface';
ss_ascii(12).txt='SAM25 R - 25 channels reverse headbox, BQ CARD Interface';
ss_ascii(13).txt='SAM32 - 32 channels headbox, Internal Interface';
ss_ascii(14).txt='SAM25 - 25 channels headbox, Internal Interface';
ss_ascii(15).txt='SAM32 R - 32 channels reverse headbox, Internal Interface';
ss_ascii(16).txt='SAM25 R - 25 channels reverse headbox, Internal Interface';
ss_ascii(18).txt='128'; %OD

fprintf(f_hdr,'%s\n','[Acquisition_Equipment]');
a_s=fread(f_trc,1,'short');
try
    a_string=['Equipment=' ss_ascii(a_s+1).txt ];
catch
    a_string=['Equipment=?'];       %OD
end
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');


% get the file type
ss_ascii(41).txt='Common reference, 128 channels of EEG';
ss_ascii(43).txt='Common reference, 84 channels of EEG, 44 channels of polygraphy';
ss_ascii(45).txt='Common reference, 84 channels of EEG, 4 reference signals (named MKR,MKRB,MKRC,MKRD)';
ss_ascii(47).txt='Common reference, 96 channels of EEG';
ss_ascii(49).txt='Common reference, 63 channels of EEG, 33 channels of polygraphy';
ss_ascii(51).txt='Common reference, 63 channels of EEG, 3 reference signals (named MKR,MKRB,MKRC)';
ss_ascii(53).txt='Common reference, 64 channels of EEG';
ss_ascii(55).txt='Common reference, 42 channels of EEG, 22 channels of polygraphy';
ss_ascii(57).txt='Common reference, 42 channels of EEG, 2 reference signals (named MKR,MKRB)';
ss_ascii(59).txt='Common reference, 32 channels of EEG';
ss_ascii(61).txt='Common reference, 21 channels of EEG, 11 channels of polygraphy';
ss_ascii(63).txt='Common reference, 21 channels of EEG, 1 reference signal (named MKR)';
ss_ascii(65).txt='Common reference, 19 channels of EEG, variable channels of polygraphy';
ss_ascii(67).txt='Common reference, 19 channels of EEG, 1 reference signal (named MKR)';
ss_ascii(69).txt='Common reference, 12 channels of EEG';
ss_ascii(71).txt='Common reference, 8 channels of EEG, variable channels of polygraphy';
ss_ascii(73).txt='Common reference, 8 channels of EEG';
ss_ascii(75).txt='Common reference, variable channels of EEG, variable channels of polygraphy';

fprintf(f_hdr,'%s\n','[File_Type]');
a_s=fread(f_trc,1,'ushort');
a_string=['Equipment=' ss_ascii(a_s+1).txt ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');


fff_raw.file_format.read_function='r_micmed';

% get adress of data
fprintf(f_hdr,'%s\n','[Adress_of_data]');
a_sad=fread(f_trc,1,'ulong');
a_string=['Adress_of_data=' num2str(a_sad) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get number of stored channels
fprintf(f_hdr,'%s\n','[Number_of_stored_channels]');
a_s=fread(f_trc,1,'ushort');
a_string=['Number_of_channels=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');
s_nchannels=a_s;

fff_raw.file_format.number_of_channels=s_nchannels;

% get distance in bytes between two successive samples
fprintf(f_hdr,'%s\n','[Distance_between_successive_samples]');
a_s=fread(f_trc,1,'ushort');
a_string=['Distance_in_bytes=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get minimum sampling rate
fprintf(f_hdr,'%s\n','[Minimum_sampling_rate]');
a_s=fread(f_trc,1,'ushort');
a_string=['Minimum_sampling_rate=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.file_format.sampling_rate=a_s;

% get size of one data
fprintf(f_hdr,'%s\n','[Size_of_each_data_in_Bytes]');
a_s=fread(f_trc,1,'ushort');
a_string=['Data_size=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.file_format.size_of_data_bytes=a_s;
s_adress_of_data=a_sad;

% get compression status
fprintf(f_hdr,'%s\n','[Compression_status]');
a_s=fread(f_trc,1,'ushort');
a_string=['Compression_status=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get number of montages used
fprintf(f_hdr,'%s\n','[Number_of_montages]');
a_s=fread(f_trc,1,'ushort');
a_string=['Number_of_montages=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get digital video start sample
fprintf(f_hdr,'%s\n','[Digital_video_start_sample]');
a_s=fread(f_trc,1,'ulong');
a_string=['Digital_video_start_sample=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

a_s=fread(f_trc,17,'uchar');

% get header type
ss_ascii(1).txt2='Micromed "System 1" Header type';
ss_ascii(2).txt2='Micromed "System 1" Header type';
ss_ascii(3).txt2='Micromed "System 2" Header Type';
ss_ascii(4).txt2='Micromed "System98" Header Type';

fprintf(f_hdr,'%s\n','[Header_Type]');
a_s=fread(f_trc,1,'uchar');
a_string=['Header_Type=' ss_ascii(a_s+1).txt2 ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');
						
% get location of info for bunch of stuff
a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.code=a_ascii;
ss_start_offset.code=fread(f_trc,1,'ulong');
ss_length.code=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.elec=a_ascii;
ss_start_offset.elec=fread(f_trc,1,'ulong');
ss_length.elec=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.note=a_ascii;
ss_start_offset.note=fread(f_trc,1,'ulong');
ss_length.note=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.flag=a_ascii;
ss_start_offset.flag=fread(f_trc,1,'ulong');
ss_length.flag=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.redu=a_ascii;
ss_start_offset.redu=fread(f_trc,1,'ulong');
ss_length.redu=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.begi=a_ascii;
ss_start_offset.begi=fread(f_trc,1,'ulong');
ss_length.begi=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.endi=a_ascii;
ss_start_offset.endi=fread(f_trc,1,'ulong');
ss_length.endi=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.mont=a_ascii;
ss_start_offset.mont=fread(f_trc,1,'ulong');
ss_length.mont=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.comp=a_ascii;
ss_start_offset.comp=fread(f_trc,1,'ulong');
ss_length.comp=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.res=a_ascii;
ss_start_offset.res=fread(f_trc,1,'ulong');
ss_length.res=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.hist=a_ascii;
ss_start_offset.hist=fread(f_trc,1,'ulong');
ss_length.hist=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.res2=a_ascii;
ss_start_offset.res2=fread(f_trc,1,'ulong');
ss_length.res2=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.eva=a_ascii;
ss_start_offset.eva=fread(f_trc,1,'ulong');
ss_length.eva=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii'
ss_name.evb=a_ascii;
ss_start_offset.evb=fread(f_trc,1,'ulong');
ss_length.evb=fread(f_trc,1,'ulong');

% get the order of the electrodes
fseek(f_trc,ss_start_offset.code,-1);
for s_c=1:s_nchannels
% 	ss_elec(s_c).order=fread(f_trc,1,'char');   
	ss_elec(s_c).order=fread(f_trc,1,'uint8');   
end; % for s_c

% get the info for each electrodes
ss_ascii(1).txt3='Referred to G2';
ss_ascii(2).txt3='Bipolar';

ss_ascii(1).txt4='nanoVolts';
ss_ascii(2).txt4='microVolts'; 
ss_ascii(3).txt4='milliVolts';
ss_ascii(4).txt4='Volts';	
ss_ascii(102).txt4= 'percents';
ss_ascii(103).txt4='bpm';
ss_ascii(104).txt4='Adim';

fseek(f_trc,ss_start_offset.elec,-1);
fprintf(f_hdr,'%s\n','');

% ss_channel will be used to store info about the channels
ss_channel.id='channel';
ss_channel.unit='none';
ss_channel.type='catego';
ss_channel.filter=ones(1,s_nchannels);

% ss_h is used to store info about each channel
for s_c=1:s_nchannels
   clear ss_h;
   
   fseek(f_trc,ss_start_offset.elec+128*ss_elec(s_c).order,-1);
   
   fprintf(f_hdr,'%s\n',['[Channel_' num2str(s_c) ']']);
   a_string=['Order=' num2str(ss_elec(s_c).order)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).type=ss_ascii(s_h+1).txt3;   
	a_string=['Type=' ss_elec(s_c).type ];
	fprintf(f_hdr,'%s\n',a_string);
   
   a_h=fread(f_trc,6,'char');
   ss_elec(s_c).p_i_l=char(a_h');
	a_string=['Positive_Input_Label=' ss_elec(s_c).p_i_l ];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.value=char(a_h');
   ss_h.authorized=1;
   ss_h.description.id=char(a_h');
   ss_h.description.positive_input_label=char(a_h');
   
   
   
   a_h=fread(f_trc,6,'char');
   ss_elec(s_c).n_i_l=char(a_h');
	a_string=['Negative_Input_Label=' ss_elec(s_c).n_i_l ];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.negative_input_label=char(a_h');
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).l_mi=s_h;   
   a_string=['Logic_Minimum=' num2str(ss_elec(s_c).l_mi)];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.logic_minimum=s_h;
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).l_ma=s_h;   
   a_string=['Logic_Maximum=' num2str(ss_elec(s_c).l_ma)];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.logic_maximum=s_h;
   
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).l_g=s_h;   
	a_string=['Logic_Ground=' num2str(ss_elec(s_c).l_g)];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.logic_ground=s_h;
   
   
   s_h=fread(f_trc,1,'long');
	ss_elec(s_c).phy_mi=s_h;   
	a_string=['Physic_Minimum=' num2str(ss_elec(s_c).phy_mi)];
   fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.physic_minimum=s_h;

   
   s_h=fread(f_trc,1,'long');
	ss_elec(s_c).phy_ma=s_h;   
	a_string=['Physic_Maximum=' num2str(ss_elec(s_c).phy_ma)];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.physic_maximum=s_h;
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).units=ss_ascii(s_h+2).txt4;   
	a_string=['Measurement_Units=' ss_elec(s_c).units ];
	fprintf(f_hdr,'%s\n',a_string);
   
   ss_h.description.measurement_units=s_h;   
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).p_h_l=s_h;   
	a_string=['Prefiltering_High_Pass_Limit=' num2str(ss_elec(s_c).p_h_l)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).p_h_t=s_h;   
	a_string=['Prefiltering_High_Type=' num2str(ss_elec(s_c).p_h_t)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).p_l_l=s_h;   
	a_string=['Prefiltering_Low_Pass_Limit=' num2str(ss_elec(s_c).p_l_l)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).p_l_t=s_h;   
	a_string=['Prefiltering_Low_Pass_Type=' num2str(ss_elec(s_c).p_l_t)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
	ss_elec(s_c).rate_coeff=s_h;   
	a_string=['Rate_Coefficient=' num2str(ss_elec(s_c).rate_coeff)];
	fprintf(f_hdr,'%s\n',a_string);
   
   s_h=fread(f_trc,1,'ushort');
   s_h=fread(f_trc,22,'char');
   a_h=fread(f_trc,64,'char');
   ss_elec(s_c).descrip=char(a_h');
	a_string=['Description=' ss_elec(s_c).descrip ];
	fprintf(f_hdr,'%s\n',a_string);
	fprintf(f_hdr,'%s\n','');

	ss_channel.values(s_c)=ss_h;

end; % for s_c

% get the notes
fseek(f_trc,ss_start_offset.note,-1);
s_h=fread(f_trc,1,'ulong');
fseek(f_trc,ss_start_offset.note,-1);
s_c=1;
clear ss_note;
while (s_h)
   s_h=fread(f_trc,1,'ulong');
	fprintf(f_hdr,'%s\n',['[Note_' num2str(s_c) ']']);
	a_string=['Time_In_Sample=' num2str(s_h)];
	fprintf(f_hdr,'%s\n',a_string);
    ss_note(s_c).time_in_sample=s_h;
    
   a_h=fread(f_trc,40,'char');
	a_string=['Text=' char(a_h')]
	fprintf(f_hdr,'%s\n',a_string);
	fprintf(f_hdr,'%s\n','');
    ss_note(s_c).text=char(a_h');
	s_c=s_c+1;
end; % while (s_h)

% get the flags
fseek(f_trc,ss_start_offset.flag,-1);
s_h=fread(f_trc,1,'ulong');
fseek(f_trc,ss_start_offset.flag,-1);
s_c=1;
while (s_h)
   s_h=fread(f_trc,1,'ulong');
   if (s_h)
      fprintf(f_hdr,'%s\n',['[Flag_' num2str(s_c) ']']);
		a_string=['Starting_Sample=' num2str(s_h)];
		fprintf(f_hdr,'%s\n',a_string);
   
   	s_h2=fread(f_trc,1,'ulong');
		a_string=['Ending_Sample=' num2str(s_h2)];
		fprintf(f_hdr,'%s\n',a_string);
		fprintf(f_hdr,'%s\n','');
      
      ss_flag(s_c).start=s_h;
      ss_flag(s_c).stop=s_h2;
      
      s_c=s_c+1;
   end; % if (s_h)   
end; % while (s_h)

% % get the segments
% fseek(f_trc,ss_start_offset.redu,-1);
% s_h=fread(f_trc,1,'ulong');
% fseek(f_trc,ss_start_offset.redu,-1);
% s_c=1;
% if (s_h==0)
%     ss_note=[];
% else
%     clear ss_note;
%     while (s_h)
%         s_h=fread(f_trc,1,'ulong');
%         fprintf(f_hdr,'%s\n',['[Note_' num2str(s_c) ']']);
%         a_string=['Time_In_Sample=' num2str(s_h)];
%         fprintf(f_hdr,'%s\n',a_string);
%         ss_note(s_c).time_in_sample=s_h;
%         
%         a_h=fread(f_trc,40,'char');
%         a_string=['Text=' char(a_h')];
%         fprintf(f_hdr,'%s\n',a_string);
%         fprintf(f_hdr,'%s\n','');
%         ss_note(s_c).text=char(a_h');
%         s_c=s_c+1;
%     end; % while (s_h)
% end;

fclose(f_trc);
fclose(f_hdr);

[v_neworder, v_neworder_rev]=mm_creanat_sys3_mono(ss_channel,a_path);
m_bipole=[];

% [v_neworder, v_neworder_rev, m_bipole]=mm_creanat_sys3(ss_channel,a_path);

% here at the very end, I have to extract the bipole matrix.
function [v_neworder, v_neworder_rev]=mm_creanat_sys3_mono(ss_channel,a_path)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
EEGFlag=1;
for s_c=1:length(ss_channel.values)
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    a_string=a_string(~(a_string=='.'));
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
%     % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
%     if (isempty(v_num))
%         s_test=0;
%     else
%         s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
%         s_test=s_test*(max(v_num)==length(a_string));
%         s_test=s_test*((min(v_num)==2)|(min(v_num)==3));
%         if (min(v_num)==2)
%             s_test=s_test*(isletter(a_string(1)));   
%             a_rt=a_string(1);
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==3)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)=='''');
%             a_rt=a_string(1:2);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;
%         
%     end;
    
    
    
        %assumes all SEEG channels are the first in the file, their last
    %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter 
    if (isempty(v_num))
        s_test=0;
%         EEGFlag=0;
    elseif (length(v_num)>2)
        s_test=0;
%         EEGFlag=0;
    elseif (max(diff(v_num))>1)
        s_test=0;
%         EEGFlag=0;
    elseif ~strcmp(deblank(ss_channel.values(s_c).description.negative_input_label),'G2')
        s_test=0;
    elseif ~EEGFlag
        s_test=0;
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
        s_test=s_test*(max(v_num)==length(a_string));
        s_test=s_test*(isletter(a_string(1)));   
        a_rt=a_string(1:(min(v_num)-1));
        s_id=str2num(a_string(v_num));
        a_name=a_string;
        if strcmp(lower(a_rt),'ecg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'ekg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'myo')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'oc')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'el')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'eog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myod')
            s_test=0;
%             EEGFlag=0;
        end
        if s_id>40
            s_test=0;
        end
            
    end

    
    
    
    
    if (s_test)
        v_find=strmatch(a_rt,v_rt,'exact');   
        if (isempty(v_find))
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
            ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)     
        end; % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
        s_other=s_other+1;   
    end; % if (s_test)   
    
end; % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end;
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    s_total=s_total+length(v_id);
end; % for s_c   

% % at this point, the electrodes that have not beed reordered should have a
% % zero in s_c
% v_find_zero=find(v_c==0);
% for s_c=1:length(v_find_zero)
%     s_ii=v_find_zero(s_c);
%     s_total=s_total+1;
%     v_c(s_ii)=s_total;
% end; % for s_c    

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);

v_neworder=v_neworder(find(v_neworder>0));

% well. it seems to work. try on other examples, allow to name the output file,
% and write a file that contains the bipoles ...


% here at the very end, I have to extract the bipole matrix.
function [v_neworder, v_neworder_rev, m_bipole]=mm_creanat_sys3(ss_channel,a_path)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
EEGFlag=1;
for s_c=1:length(ss_channel.values)
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    a_string=a_string(~(a_string=='.'));
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
%     % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
%     if (isempty(v_num))
%         s_test=0;
%     else
%         s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
%         s_test=s_test*(max(v_num)==length(a_string));
%         s_test=s_test*((min(v_num)==2)|(min(v_num)==3));
%         if (min(v_num)==2)
%             s_test=s_test*(isletter(a_string(1)));   
%             a_rt=a_string(1);
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==3)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)=='''');
%             a_rt=a_string(1:2);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;
%         
%     end;
    
    
    
        %assumes all SEEG channels are the first in the file, their last
    %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter 
    if (isempty(v_num))
        s_test=0;
%         EEGFlag=0;
    elseif (length(v_num)>2)
        s_test=0;
%         EEGFlag=0;
    elseif (max(diff(v_num))>1)
        s_test=0;
%         EEGFlag=0;
    elseif ~strcmp(deblank(ss_channel.values(s_c).description.negative_input_label),'G2')
        s_test=0;
    elseif ~EEGFlag
        s_test=0;
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
        s_test=s_test*(max(v_num)==length(a_string));
        s_test=s_test*(isletter(a_string(1)));   
        a_rt=a_string(1:(min(v_num)-1));
        s_id=str2num(a_string(v_num));
        a_name=a_string;
        if strcmp(lower(a_rt),'ecg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'ekg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'myo')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'oc')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'el')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'eog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myod')
            s_test=0;
%             EEGFlag=0;
        end
        if s_id>40
            s_test=0;
        end
            
    end

    
    
    
    
    if (s_test)
        v_find=strmatch(a_rt,v_rt,'exact');   
        if (isempty(v_find))
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
            ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)     
        end; % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
        s_other=s_other+1;   
    end; % if (s_test)   
    
end; % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
m_bipole=[];
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end;
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    % we also want to create the bipole txt file
    % that means we want within each electrode, find the site pairs of the type i i+1
    for s_b=2:length(v_y)
        if ((v_y(s_b)-v_y(s_b-1))==1) % then we have a bipole
            m_bipole=[m_bipole;[s_total+s_b s_total+s_b-1]];         
        end;
    end;   
    s_total=s_total+length(v_id);
end; % for s_c   

% at this point, the electrodes that have not beed reordered should have a
% zero in s_c
v_find_zero=find(v_c==0);
for s_c=1:length(v_find_zero)
    s_ii=v_find_zero(s_c);
    s_total=s_total+1;
    v_c(s_ii)=s_total;
end; % for s_c    

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);

f=fopen([a_path filesep 'bipole.txt'],'w');
for s_b=1:size(m_bipole,1)
    fprintf(f,'%s\t%s\n',ss_channel.values(v_neworder_rev(m_bipole(s_b,1))).value,ss_channel.values(v_neworder_rev(m_bipole(s_b,2))).value);
end;
fclose(f);

% well. it seems to work. try on other examples, allow to name the output file,
% and write a file that contains the bipoles ...

% the bipole thing should be easy inside a v_rt with the v_id ....
% then insert into lyonstim, to make lyonstim2, then convert into C
% and test it to see if it works.
% then convert ella into C ?!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fff_raw, ss_channel, s_adress_of_data, v_neworder, v_neworder_rev, m_bipole, ss_note]=rdhd_sys4(a_trcfilename, FileOut)

% This function reads the header of a micromed .trc file provided that this
% header is of type 4
% It creates a corresponding hdr file, (ascii) with all this info

% exple: [fff_raw, ss_channel, s_adress_of_data,m_bipole,ss_note]=rdhd_sys4('d:/data2005/eeg_40.trc');

% returns a structure fff_raw that has the characteristics of a fff_raw ella data to describe the exam
% in this directory, we have created a directory a_dirname with the name of the patient and the data of recording
% that's where this hdr file will go.

% get file name and path

[a_path,a_name,a_ext]=fileparts(a_trcfilename);
a_dirname=a_path;

% get the file size
a_com=['ss_d=dir(''' a_trcfilename ''');'];
eval(a_com);
file_size=ss_d.bytes;

% open file
f_trc=fopen(a_trcfilename,'r','n','windows-1252');
a_trcfilename(end-2:end)=lower(a_trcfilename(end-2:end));
a_hdrname=[FileOut '.hdr'];
f_hdr=fopen(a_hdrname,'w');

fprintf(f_hdr,'%s\n','[File_Size]');
a_string=['File_Size_in_Bytes=' num2str(file_size)];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get title
fprintf(f_hdr,'%s\n','[Title]');
a_s=fread(f_trc,32,'char');
a_string=['Title=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get lab
fprintf(f_hdr,'%s\n','[Lab]');
a_s=fread(f_trc,32,'char');
a_string=['Lab=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get patient info
fprintf(f_hdr,'%s\n','[Patient]');
a_s=fread(f_trc,22,'char');
a_string=['Last_Name=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.name=deblank(char(a_s'));

a_s=fread(f_trc,20,'char');
a_string=['First_Name=' char(a_s') ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.surname=deblank(char(a_s'));

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Month=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Day=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Birth_Year=' num2str(1900+a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fff_raw.subject.yob=(1900+a_s);

a_s=fread(f_trc,19,'char');
fprintf(f_hdr,'%s\n','');

fff_raw.subject.handedness='?'; % not in the header
fff_raw.subject.sex='?'; % not in the header

% recording date
fprintf(f_hdr,'%s\n','[Recording Date and Time]');

a_s2=fread(f_trc,1,'uchar');
a_string=['Recording_Day=' num2str(a_s2) ];
fprintf(f_hdr,'%s\n',a_string);

a_s1=fread(f_trc,1,'uchar');
a_string=['Recording_Month=' num2str(a_s1) ];
fprintf(f_hdr,'%s\n',a_string);

a_s3=fread(f_trc,1,'uchar');
a_string=['Recording_Year=' num2str(1900+a_s3) ];
fprintf(f_hdr,'%s\n',a_string);

fff_raw.recording_date.year=1900+a_s3;
fff_raw.recording_date.month=a_s1;
fff_raw.recording_date.day=a_s2;

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Hour=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Min=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);

a_s=fread(f_trc,1,'uchar');
a_string=['Recording_Sec=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.recording_site.lab='Grenoble';
fff_raw.recording_site.pi='JP_Lachaux';
fff_raw.exam=['exp_' fff_raw.subject.name '_'  int2str(fff_raw.recording_date.day) '_' int2str(fff_raw.recording_date.month) '_' int2str(fff_raw.recording_date.year)];


% get the acquisition equipment
ss_ascii(1).txt='BQ124 - 24 channels headbox, Internal Interface';
ss_ascii(2).txt='MS40 - Holter recorder';
ss_ascii(3).txt='BQ132S - 32 channels headbox, Internal Interface';
ss_ascii(7).txt='BQ124 - 24 channels headbox, BQ CARD Interface';
ss_ascii(8).txt='SAM32 - 32 channels headbox, BQ CARD Interface';
ss_ascii(9).txt='SAM25 - 25 channels headbox, BQ CARD Interface';
ss_ascii(10).txt='BQ132S R - 32 channels reverse headbox, Internal Interface';
ss_ascii(11).txt='SAM32 R - 32 channels reverse headbox, BQ CARD Interface';
ss_ascii(12).txt='SAM25 R - 25 channels reverse headbox, BQ CARD Interface';
ss_ascii(13).txt='SAM32 - 32 channels headbox, Internal Interface';
ss_ascii(14).txt='SAM25 - 25 channels headbox, Internal Interface';
ss_ascii(15).txt='SAM32 R - 32 channels reverse headbox, Internal Interface';
ss_ascii(16).txt='SAM25 R';
ss_ascii(17).txt='SD - 32 channels headbox with jackbox ';
ss_ascii(18).txt='SD128';
ss_ascii(19).txt='SD96';
ss_ascii(20).txt='SD64';
ss_ascii(21).txt='SD128c';
ss_ascii(22).txt='SD64c';
ss_ascii(23).txt='BQ132S';
ss_ascii(23).txt='BQ132SR';
ss_ascii(24).txt='BQ132SR';

ss_ascii(33).txt='SDMRI32';
% get the file type
ss_ascii(41).txt='Common reference, 128 channels of EEG';
ss_ascii(43).txt='Common reference, 84 channels of EEG, 44 channels of polygraphy';
ss_ascii(45).txt='Common reference, 84 channels of EEG, 4 reference signals (named MKR,MKRB,MKRC,MKRD)';
ss_ascii(47).txt='Common reference, 96 channels of EEG';
ss_ascii(49).txt='Common reference, 63 channels of EEG, 33 channels of polygraphy';
ss_ascii(51).txt='Common reference, 63 channels of EEG, 3 reference signals (named MKR,MKRB,MKRC)';
ss_ascii(53).txt='Common reference, 64 channels of EEG';
ss_ascii(55).txt='Common reference, 42 channels of EEG, 22 channels of polygraphy';
ss_ascii(57).txt='Common reference, 42 channels of EEG, 2 reference signals (named MKR,MKRB)';
ss_ascii(59).txt='Common reference, 32 channels of EEG';
ss_ascii(61).txt='Common reference, 21 channels of EEG, 11 channels of polygraphy';
ss_ascii(63).txt='Common reference, 21 channels of EEG, 1 reference signal (named MKR)';
ss_ascii(65).txt='Common reference, 19 channels of EEG, variable channels of polygraphy';
ss_ascii(67).txt='Common reference, 19 channels of EEG, 1 reference signal (named MKR)';
ss_ascii(69).txt='Common reference, 12 channels of EEG';
ss_ascii(71).txt='Common reference, 8 channels of EEG, variable channels of polygraphy';
ss_ascii(73).txt='Common reference, 8 channels of EEG';
ss_ascii(75).txt='Common reference, variable channels of EEG, variable channels of polygraphy';
ss_ascii(77).txt='JPL : not used';
ss_ascii(79).txt='JPL : not used';
ss_ascii(81).txt='JPL : not used';
ss_ascii(83).txt='JPL : not used';
ss_ascii(85).txt='JPL : not used';
ss_ascii(87).txt='JPL : not used';
ss_ascii(101).txt='JPL : not used';
ss_ascii(102).txt='JPL : not used';
ss_ascii(103).txt='JPL : not used';
ss_ascii(104).txt='JPL : not used';
ss_ascii(121).txt='JPL : not used';
ss_ascii(122).txt='JPL : not used';
ss_ascii(123).txt='JPL : not used';
ss_ascii(141).txt='JPL : not used';
ss_ascii(142).txt='JPL : not used';
ss_ascii(161).txt='JPL : not used';
ss_ascii(162).txt='JPL : not used';
ss_ascii(163).txt='JPL : not used';
ss_ascii(181).txt='JPL : not used';
ss_ascii(182).txt='JPL : not used';
ss_ascii(183).txt='JPL : not used';
ss_ascii(184).txt='JPL : not used';
ss_ascii(201).txt='JPL : not used';
ss_ascii(202).txt='JPL : not used';
ss_ascii(203).txt='JPL : not used';
ss_ascii(204).txt='JPL : not used';
ss_ascii(205).txt='JPL : not used';
ss_ascii(206).txt='JPL : not used';

fprintf(f_hdr,'%s\n','[Acquisition_Equipment]');
a_s=fread(f_trc,1,'short');
%a_s=0; %%% JP debug
a_string=['Equipment=' ss_ascii(a_s+1).txt ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');


fprintf(f_hdr,'%s\n','[File_Type]');
a_s=fread(f_trc,1,'ushort');
a_string=['Equipment=' ss_ascii(a_s+1).txt ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');


fff_raw.file_format.read_function='r_micmed';

% get adress of data
fprintf(f_hdr,'%s\n','[Adress_of_data]');
a_sad=fread(f_trc,1,'ulong');
a_string=['Adress_of_data=' num2str(a_sad) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get number of stored channels
fprintf(f_hdr,'%s\n','[Number_of_stored_channels]');
a_s=fread(f_trc,1,'ushort');
a_string=['Number_of_channels=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');
s_nchannels=a_s;

fff_raw.file_format.number_of_channels=s_nchannels;

% get distance in bytes between two successive samples
fprintf(f_hdr,'%s\n','[Distance_between_successive_samples]');
a_s=fread(f_trc,1,'ushort');
a_string=['Distance_in_bytes=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get minimum sampling rate
fprintf(f_hdr,'%s\n','[Minimum_sampling_rate]');
a_s=fread(f_trc,1,'ushort');
a_string=['Minimum_sampling_rate=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.file_format.sampling_rate=a_s;

% get size of one data
fprintf(f_hdr,'%s\n','[Size_of_each_data_in_Bytes]');
a_s=fread(f_trc,1,'ushort');
a_string=['Data_size=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

fff_raw.file_format.size_of_data_bytes=a_s;
s_adress_of_data=a_sad;

% get compression status
fprintf(f_hdr,'%s\n','[Compression_status]');
a_s=fread(f_trc,1,'ushort');
a_string=['Compression_status=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get number of montages used
fprintf(f_hdr,'%s\n','[Number_of_montages]');
a_s=fread(f_trc,1,'ushort');
a_string=['Number_of_montages=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get digital video start sample
fprintf(f_hdr,'%s\n','[Digital_video_start_sample]');
a_s=fread(f_trc,1,'ulong');
a_string=['Digital_video_start_sample=' num2str(a_s) ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

a_s=fread(f_trc,17,'uchar');

% get header type
ss_ascii(1).txt2='Micromed "System 1" Header type';
ss_ascii(2).txt2='Micromed "System 1" Header type';
ss_ascii(3).txt2='Micromed "System 2" Header Type';
ss_ascii(4).txt2='Micromed "System98" Header Type';
ss_ascii(5).txt2='Micromed "System98" Header Type';

fprintf(f_hdr,'%s\n','[Header_Type]');
a_s=fread(f_trc,1,'uchar');
if (a_s~=4)
    disp('Watch Out : You are NOT reading the correct Header Type ');
end;   
a_string=['Header_Type=' ss_ascii(a_s+1).txt2 ];
fprintf(f_hdr,'%s\n',a_string);
fprintf(f_hdr,'%s\n','');

% get location of info for bunch of stuff
a_ascii=char(fread(f_trc,8,'char'));
%a_ascii'
ss_name.code=a_ascii;
ss_start_offset.code=fread(f_trc,1,'ulong');
ss_length.code=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
LABCOD=a_ascii';
ss_name.elec=a_ascii;
ss_start_offset.elec=fread(f_trc,1,'ulong');
ss_length.elec=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
NOTE=a_ascii';
ss_name.note=a_ascii;
ss_start_offset.note=fread(f_trc,1,'ulong');
ss_length.note=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
FLAGS=a_ascii';
ss_name.flag=a_ascii;
ss_start_offset.flag=fread(f_trc,1,'ulong');
ss_length.flag=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
TRONCA=a_ascii';
ss_name.redu=a_ascii;
ss_start_offset.redu=fread(f_trc,1,'ulong');
ss_length.redu=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.begi=a_ascii;
ss_start_offset.begi=fread(f_trc,1,'ulong');
ss_length.begi=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.endi=a_ascii;
ss_start_offset.endi=fread(f_trc,1,'ulong');
ss_length.endi=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.mont=a_ascii;
ss_start_offset.mont=fread(f_trc,1,'ulong');
ss_length.mont=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.comp=a_ascii;
ss_start_offset.comp=fread(f_trc,1,'ulong');
ss_length.comp=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.res=a_ascii;
ss_start_offset.res=fread(f_trc,1,'ulong');
ss_length.res=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.hist=a_ascii;
ss_start_offset.hist=fread(f_trc,1,'ulong');
ss_length.hist=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.res2=a_ascii;
ss_start_offset.res2=fread(f_trc,1,'ulong');
ss_length.res2=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.eva=a_ascii;
ss_start_offset.eva=fread(f_trc,1,'ulong');
ss_length.eva=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.evb=a_ascii;
ss_start_offset.evb=fread(f_trc,1,'ulong');
ss_length.evb=fread(f_trc,1,'ulong');

a_ascii=char(fread(f_trc,8,'char'));
a_ascii';
ss_name.trig=a_ascii;
ss_start_offset.trig=fread(f_trc,1,'ulong');
ss_length.trig=fread(f_trc,1,'ulong');

% get the order of the electrodes
fseek(f_trc,ss_start_offset.code,-1);
for s_c=1:s_nchannels
    %ss_elec(s_c).order=fread(f_trc,1,'char');   %%% CHECK THAT OUT : replace with 'ushort' ?
    ss_elec(s_c).order=fread(f_trc,1,'ushort');   %%% CHECK THAT OUT : replace with 'ushort' ?
end; % for s_c

% get the info for each electrodes
ss_ascii(1).txt3='Referred to G2';
ss_ascii(2).txt3='Bipolar';

ss_ascii(1).txt4='nanoVolts';
ss_ascii(2).txt4='microVolts'; 
ss_ascii(3).txt4='milliVolts';
ss_ascii(4).txt4='Volts';	
ss_ascii(102).txt4= 'percents';
ss_ascii(103).txt4='bpm';
ss_ascii(104).txt4='Adim';

fseek(f_trc,ss_start_offset.elec,-1);
fprintf(f_hdr,'%s\n','');

% ss_channel will be used to store info about the channels
ss_channel.id='channel';
ss_channel.unit='none';
ss_channel.type='catego';
ss_channel.filter=ones(1,s_nchannels);

% ss_h is used to store info about each channel
for s_c=1:s_nchannels
    clear ss_h;
    
    fseek(f_trc,ss_start_offset.elec+128*ss_elec(s_c).order,-1);
    
    fprintf(f_hdr,'%s\n',['[Channel_' num2str(s_c) ']']);
    a_string=['Order=' num2str(ss_elec(s_c).order)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'uchar'); % add another fread for status ?
    
    s_h=fread(f_trc,1,'uchar'); % change into uchar ? add another fread for status ?
    ss_elec(s_c).type=ss_ascii(s_h+1).txt3;   
    a_string=['Type=' ss_elec(s_c).type ];
    fprintf(f_hdr,'%s\n',a_string);
    
    a_h=fread(f_trc,6,'char');
    ss_elec(s_c).p_i_l=char(a_h');
    a_string=['Positive_Input_Label=' ss_elec(s_c).p_i_l ];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.value=char(a_h');
    ss_h.authorized=1;
    ss_h.description.id=char(a_h');
    ss_h.description.positive_input_label=char(a_h');
    
    
    
    a_h=fread(f_trc,6,'char');
    ss_elec(s_c).n_i_l=char(a_h');
    a_string=['Negative_Input_Label=' ss_elec(s_c).n_i_l ];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.negative_input_label=char(a_h');
    
    s_h=fread(f_trc,1,'long');
    ss_elec(s_c).l_mi=s_h;   
    a_string=['Logic_Minimum=' num2str(ss_elec(s_c).l_mi)];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.logic_minimum=s_h;
    
    s_h=fread(f_trc,1,'long');
    ss_elec(s_c).l_ma=s_h;   
    a_string=['Logic_Maximum=' num2str(ss_elec(s_c).l_ma)];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.logic_maximum=s_h;
    
    
    s_h=fread(f_trc,1,'long');
    ss_elec(s_c).l_g=s_h;   
    a_string=['Logic_Ground=' num2str(ss_elec(s_c).l_g)];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.logic_ground=s_h;
    
    
    s_h=fread(f_trc,1,'long');
    ss_elec(s_c).phy_mi=s_h;   
    a_string=['Physic_Minimum=' num2str(ss_elec(s_c).phy_mi)];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.physic_minimum=s_h;
    
    
    s_h=fread(f_trc,1,'long');
    ss_elec(s_c).phy_ma=s_h;   
    a_string=['Physic_Maximum=' num2str(ss_elec(s_c).phy_ma)];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.physic_maximum=s_h;
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).units=ss_ascii(s_h+2).txt4;   
    a_string=['Measurement_Units=' ss_elec(s_c).units ];
    fprintf(f_hdr,'%s\n',a_string);
    
    ss_h.description.measurement_units=s_h;   
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).p_h_l=s_h;   
    a_string=['Prefiltering_High_Pass_Limit=' num2str(ss_elec(s_c).p_h_l)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).p_h_t=s_h;   
    a_string=['Prefiltering_High_Type=' num2str(ss_elec(s_c).p_h_t)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).p_l_l=s_h;   
    a_string=['Prefiltering_Low_Pass_Limit=' num2str(ss_elec(s_c).p_l_l)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).p_l_t=s_h;   
    a_string=['Prefiltering_Low_Pass_Type=' num2str(ss_elec(s_c).p_l_t)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'ushort');
    ss_elec(s_c).rate_coeff=s_h;   
    a_string=['Rate_Coefficient=' num2str(ss_elec(s_c).rate_coeff)];
    fprintf(f_hdr,'%s\n',a_string);
    
    s_h=fread(f_trc,1,'ushort');
    s_h=fread(f_trc,1,'float');
    s_h=fread(f_trc,1,'float');
    s_h=fread(f_trc,1,'uchar');
    s_h=fread(f_trc,1,'uchar');
    a_h=fread(f_trc,32,'char');
    ss_elec(s_c).descrip=char(a_h');
    a_string=['Description=' ss_elec(s_c).descrip ];
    fprintf(f_hdr,'%s\n',a_string);
    fprintf(f_hdr,'%s\n','');
    
    s_h=fread(f_trc,1,'float');
    s_h=fread(f_trc,1,'float');
    s_h=fread(f_trc,1,'float');
    
    ss_channel.values(s_c)=ss_h;
    
end; % for s_c

% get the notes
fseek(f_trc,ss_start_offset.note,-1);
s_h=fread(f_trc,1,'ulong');
fseek(f_trc,ss_start_offset.note,-1);
s_c=1;
if (s_h==0)
    ss_note=[];
else
    clear ss_note;
    while (s_h)
        s_h=fread(f_trc,1,'ulong');
        fprintf(f_hdr,'%s\n',['[Note_' num2str(s_c) ']']);
        a_string=['Time_In_Sample=' num2str(s_h)];
        fprintf(f_hdr,'%s\n',a_string);
        ss_note(s_c).time_in_sample=s_h;
        
        a_h=fread(f_trc,40,'char');
        a_string=['Text=' char(a_h')];
        fprintf(f_hdr,'%s\n',a_string);
        fprintf(f_hdr,'%s\n','');
        ss_note(s_c).text=char(a_h');
        s_c=s_c+1;
    end; % while (s_h)
end;

fclose(f_trc);
fclose(f_hdr);

[v_neworder, v_neworder_rev, ss_channel]=mm_creanat_sys4_mono(ss_channel,a_path);
m_bipole=[];

% switch Montage
%     case 'Gre Epi'
%         [v_neworder, v_neworder_rev, m_bipole, ss_channel]=mm_creanat_sys4(ss_channel,a_path);
%     case 'Gre Park'
%         [v_neworder, v_neworder_rev, m_bipole]=mm_creanat_sys4_Park(ss_channel,a_path);
% end


% here at the very end, I have to extract the bipole matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_neworder, v_neworder_rev, ss_channel]=mm_creanat_sys4_mono(ss_channel,a_path)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
EEGFlag=1;
for s_c=1:length(ss_channel.values)
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    
    %OD for Lyon
    a_string=lower(a_string(~isspace(a_string)));
    %OD for Brousse
    a_string=a_string(~(a_string=='.'));
    %OD for StA12Pet
    if strcmp('StA12Pet',spm_str_manip(a_path,'t'))
        if strcmp(a_string,'g''1')
            a_string='g1';
        elseif strcmp(a_string,'g''2')
            a_string='g2';
        end
    end

    ss_channel.values(s_c).value=a_string;

    
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    
%     % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
%     if (isempty(v_num))
%         s_test=0;
%     else
%         s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
%         s_test=s_test*(max(v_num)==length(a_string));
%         s_test=s_test*((min(v_num)==2)|(min(v_num)==3)|(min(v_num)==4));
%         if (min(v_num)==2)
%             s_test=s_test*(isletter(a_string(1)));   
%             a_rt=a_string(1);
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==3)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)==''''||a_string(2)=='t');
%             a_rt=a_string(1:2);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==4)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)=='t');
%             s_test=s_test*(a_string(3)=='''');
%             a_rt=a_string(1:3);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%     end; % if isempty
    
    
    %assumes all SEEG channels are the first in the file, their last
    %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter 
    if (isempty(v_num))
        s_test=0;
%         EEGFlag=0;
    elseif (length(v_num)>2)
        s_test=0;
%         EEGFlag=0;
    elseif (max(diff(v_num))>1)
        s_test=0;
%         EEGFlag=0;
    elseif ~strcmp(deblank(ss_channel.values(s_c).description.negative_input_label),'G2')
        s_test=0;
    elseif ~EEGFlag
        s_test=0;
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
        s_test=s_test*(max(v_num)==length(a_string));
        s_test=s_test*(isletter(a_string(1)));   
        a_rt=a_string(1:(min(v_num)-1));
        s_id=str2num(a_string(v_num));
        a_name=a_string;
        if strcmp(lower(a_rt),'ecg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'ekg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'myo')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'oc')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'dd')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'dg')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'el')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'eog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myod')
            s_test=0;
%             EEGFlag=0;
        end
        if s_id>40
            s_test=0;
        end
            
    end
        
    if (s_test)
        v_find=strmatch(a_rt,v_rt,'exact');   
        if (isempty(v_find))
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            %avoid the possibility to have two electrodes with the same name (happens in Nancy): Take the first occurence
            if isempty(intersect(ss_elec(s_count).v_id,s_id))
                ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
                ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)
            end
        end; % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
        s_other=s_other+1;   
    end; % if (s_test)   
end; % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end;
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    s_total=s_total+length(v_id);
end; % for s_c   

% % at this point, the electrodes that have not beed reordered should have a
% % zero in s_c
% v_find_zero=find(v_c==0);
% for s_c=1:length(v_find_zero)
%     s_ii=v_find_zero(s_c);
%     s_total=s_total+1;
%     v_c(s_ii)=s_total;
% end; % for s_c    

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);

v_neworder=v_neworder(find(v_neworder>0));

% well. it seems to work. try on other examples, allow to name the output file,
% and write a file that contains the bipoles ...

% the bipole thing should be easy inside a v_rt with the v_id ....
% then insert into lyonstim, to make lyonstim2, then convert into C
% and test it to see if it works.
% then convert ella into C ?!



% here at the very end, I have to extract the bipole matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_neworder, v_neworder_rev, m_bipole, ss_channel]=mm_creanat_sys4(ss_channel,a_path,Montage)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
EEGFlag=1;
for s_c=1:length(ss_channel.values)
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    
    %OD for Lyon
    a_string=lower(a_string(~isspace(a_string)));
    %OD for Brousse
    a_string=a_string(~(a_string=='.'));
    %OD for StA12Pet
    if strcmp('StA12Pet',spm_str_manip(a_path,'t'))
        if strcmp(a_string,'g''1')
            a_string='g1';
        elseif strcmp(a_string,'g''2')
            a_string='g2';
        end
    end

    ss_channel.values(s_c).value=a_string;

    
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    
%     % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
%     if (isempty(v_num))
%         s_test=0;
%     else
%         s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
%         s_test=s_test*(max(v_num)==length(a_string));
%         s_test=s_test*((min(v_num)==2)|(min(v_num)==3)|(min(v_num)==4));
%         if (min(v_num)==2)
%             s_test=s_test*(isletter(a_string(1)));   
%             a_rt=a_string(1);
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==3)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)==''''||a_string(2)=='t');
%             a_rt=a_string(1:2);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%         if (min(v_num)==4)
%             s_test=s_test*(isletter(a_string(1)));
%             s_test=s_test*(a_string(2)=='t');
%             s_test=s_test*(a_string(3)=='''');
%             a_rt=a_string(1:3);   
%             s_id=str2num(a_string(v_num));
%             a_name=a_string;   
%         end;   
%     end; % if isempty
    
    
    %assumes all SEEG channels are the first in the file, their last
    %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter 
    if (isempty(v_num))
        s_test=0;
%         EEGFlag=0;
    elseif (length(v_num)>2)
        s_test=0;
%         EEGFlag=0;
    elseif (max(diff(v_num))>1)
        s_test=0;
%         EEGFlag=0;
    elseif ~strcmp(deblank(ss_channel.values(s_c).description.negative_input_label),'G2')
        s_test=0;
    elseif ~EEGFlag
        s_test=0;
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
        s_test=s_test*(max(v_num)==length(a_string));
        s_test=s_test*(isletter(a_string(1)));   
        a_rt=a_string(1:(min(v_num)-1));
        s_id=str2num(a_string(v_num));
        a_name=a_string;
        if strcmp(lower(a_rt),'ecg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'ekg')
            s_test=1;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'myo')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'oc')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'dd')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'dg')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'el')
            s_test=0;
%             EEGFlag=0;
        end    
        if strcmp(lower(a_rt),'eog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myog')
            s_test=0;
%             EEGFlag=0;
        end
        if strcmp(lower(a_rt),'myod')
            s_test=0;
%             EEGFlag=0;
        end
        if s_id>40
            s_test=0;
        end
            
    end
        
    if (s_test)
        v_find=strmatch(a_rt,v_rt,'exact');   
        if (isempty(v_find))
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
            ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)     
        end; % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
        s_other=s_other+1;   
    end; % if (s_test)   
end; % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
m_bipole=[];
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end;
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    % we also want to create the bipole txt file
    % that means we want within each electrode, find the site pairs of the type i i+1
    for s_b=2:length(v_y)
        if ((v_y(s_b)-v_y(s_b-1))==1) % then we have a bipole
            m_bipole=[m_bipole;[s_total+s_b s_total+s_b-1]];         
        end;
    end;   
    s_total=s_total+length(v_id);
end; % for s_c   

% at this point, the electrodes that have not beed reordered should have a
% zero in s_c
v_find_zero=find(v_c==0);
for s_c=1:length(v_find_zero)
    s_ii=v_find_zero(s_c);
    s_total=s_total+1;
    v_c(s_ii)=s_total;
end; % for s_c    

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);

f=fopen([a_path filesep 'bipole.txt'],'w');
for s_b=1:size(m_bipole,1)
    fprintf(f,'%s\t%s\n',ss_channel.values(v_neworder_rev(m_bipole(s_b,1))).value,ss_channel.values(v_neworder_rev(m_bipole(s_b,2))).value);
end;
fclose(f);

% well. it seems to work. try on other examples, allow to name the output file,
% and write a file that contains the bipoles ...

% the bipole thing should be easy inside a v_rt with the v_id ....
% then insert into lyonstim, to make lyonstim2, then convert into C
% and test it to see if it works.
% then convert ella into C ?!


% here at the very end, I have to extract the bipole matrix.
function [v_neworder, v_neworder_rev, m_bipole]=mm_creanat_sys4_Park(ss_channel,a_path)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
for s_c=1:length(ss_channel.values)
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    a_namex=a_string;
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
    if (isempty(v_num))
        s_test=0;
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
        s_test=s_test*(max(v_num)==length(a_string));
        s_test=s_test*((min(v_num)==2)|(min(v_num)==3) | (strcmp(a_namex(1:min(3,length(a_namex))),'stn')));
        if (min(v_num)==2)
            s_test=s_test*(isletter(a_string(1)));   
            a_rt=a_string(1);
            s_id=str2num(a_string(v_num));
            a_name=a_string;   
        end;   
        if (min(v_num)==3)
            s_test=s_test*(isletter(a_string(1)));
            s_test=s_test*(a_string(2)=='''');
            a_rt=a_string(1:2);   
            s_id=str2num(a_string(v_num));
            a_name=a_string;
        end;
        if (strcmp(a_namex(1:min(3,length(a_namex))),'stn'))
            s_test=1;
            a_rt=a_string(1:3);   
            s_id=str2num(a_string(v_num));
            a_name=a_string;   
        end;   
    end; % if isempty
    
    if (s_test)
        v_find=strmatch(a_rt,v_rt,'exact');   
        if (isempty(v_find))
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
            ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)     
        end; % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)     
        s_other=s_other+1;   
    end; % if (s_test)   
end; % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
m_bipole=[];
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end;
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    % we also want to create the bipole txt file
    % that means we want within each electrode, find the site pairs of the type i i+1
    for s_b=2:length(v_y)
        if ((v_y(s_b)-v_y(s_b-1))==1) % then we have a bipole
            m_bipole=[m_bipole;[s_total+s_b s_total+s_b-1]];         
        end;
    end;   
    s_total=s_total+length(v_id);
end; % for s_c   

% at this point, the electrodes that have not beed reordered should have a
% zero in s_c
v_find_zero=find(v_c==0);
for s_c=1:length(v_find_zero)
    s_ii=v_find_zero(s_c);
    s_total=s_total+1;
    v_c(s_ii)=s_total;
end; % for s_c    

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);

%Pour Julien
%A AMELIORER!!!!!!!!!
% v_neworder_rev=[1:21]';
% m_bipole=[1 0
%     2 0
%     3 0
%     4 0
%     5 0
%     6 0
%     7 0
%     8 0
%     9 0
%     10 0
%     11 0
%     12 0
%     13 0
%     15 14
%     16 15
%     17 16
%     19 18
%     20 19
%     21 20];


v_neworder_rev=[1:17]';
m_bipole=[1 0
    2 0
    3 0
    4 0
    5 0
    6 0
    7 0
    8 0
    9 0
    11 10
    12 11
    13 12
    15 14
    16 15
    17 16];

f=fopen([pwd filesep 'bipole.txt'],'w');
for s_b=1:size(m_bipole,1)
    if m_bipole(s_b,2)~=0
        fprintf(f,'%s\t%s\n',ss_channel.values(v_neworder_rev(m_bipole(s_b,1))).value,ss_channel.values(v_neworder_rev(m_bipole(s_b,2))).value);
    else
        fprintf(f,'%s\n',ss_channel.values(v_neworder_rev(m_bipole(s_b,1))).value);
    end
end;
fclose(f);

% well. it seems to work. try on other examples, allow to name the output file,
% and write a file that contains the bipoles ...

% the bipole thing should be easy inside a v_rt with the v_id ....
% then insert into lyonstim, to make lyonstim2, then convert into C
% and test it to see if it works.
% then convert ella into C ?!










% [DATAOUT]=readtrc(SETTIINGS)  -   Reads Micromed System 98 *.trc EEG File 
%                                   for import into EEGLAB. This file reads 
%                                   Micromed System Plus EEG *.trc files with 
%                                   header of type 4
%
% USAGE:
%   >> [DATAOUT]=readtrc(SETTINGS)
%
%
% INPUT
%   SETTINGS is a struct holding the parameters for reading the .TRC file
%   SETTINGS has the following fields:
%
%       SETTINGS.filename :                 Name of file to be imported
%       SETTINGS.loadevents.state :         'yes' for loading event triggers
%                                           'no' for not
%       SETTINGS.loadevents.type :          'marker' for event triggers inserted 
%                                           on 'MKR channel
%                                           'eegchan' for triggers inserted on 
%                                           EEG channels
%                                           'none' or
%                                           'both'
%       SETTINGS.loadevents.dig_ch1:        number of name of eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch1_label:  label to give events on eegchan marker channel 1
%       SETTINGS.loadevents.dig_ch2:        number of name of eegchan marker channel 2
%       SETTINGS.loadevents.dig_ch2_label:  label to give events on eegchan marker channel 2
%       SETTINGS.chan_adjust_status:        1 for adjusting amp of channels 0 for not
%       SETTINGS.chans                      channels to load, [ ] for all
%       (default)
%       SETTINGS.chan_adjust                channels to adjust
%
%   Alternant method: enter 11 input arguments each corresponding to one of the above
%   fields.  If only one arg is used,  it must be the struct above.  If more, there
%   must be 11 inputs in the order above i.e. OUT=readtrc(filename,eventstate....etc).
%
% OUTPUT
%   DATAOUT with same fields as EEGLAB EEG file structure. (see eeg_checkset.m)
%   The following fields are use: 
%       DATAOUT.data 
%       DATAOUT.filename
%       DATAOUT.filepath
%       DATAOUT.srate
%       DATAOUT.setname
%       DATAOUT.pnts
%       DATAOUT.nbchan
%       DATAOUT.trials
%       DATAOUT.xmin
%       DATAOUT.ref
%
% Author: Rami K. Niazy
% Copyright (c) University of Oxford.
%
%   see also pop_readtrc.m    eegplugin_trcimport.m

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004-2005 University of Oxford
% Author: Rami K. Niazy
%         rami@fmrib.ox.ac.uk
%
% This program can only be provided from Micromed s.r.l.  You may not
% give away or edit this program. 

% April 19, 2005
% Version 1.1
% Reads triggers from 'MRK' and EEG channels together
% Assigns value of 'MRK' trigger as event 'type' in EEGLAB
% Reads only selected channels into memory

% Dec 23, 2004
% Fixed bug for No MRK case

% Sep 29, 2004
% Fixed Trigger offset

% Sep 22, 2004
% Fixed Trigger offset

% Sep 21, 2004
% Allow usage of ':' in specifying
% channels to exclude or adjust

% Aug 12, 2004
% Added command line usage
% fixed com out for EEGLAB

% Aug 10, 2004
% Adjusted Licence
% Added Chan exclude
% Fixed scaling

% Aug 4, 2004
% Fixed scaling of data to uV
% Fixed Analog Trig size detection

% Date: Julyl 30, 2004
% Initial Setup

function [TRC]=readtrc98(varargin)

if nargin<1
    error('Not enough input arguments. Please see help file for usage');
elseif nargin==1
    if isstruct(varargin{1})
        PARAM=varargin{1};
    else
        error('Single input usage must be a Structure.  Please see help file.');
    end
elseif nargin < 10
    error('Not enough input arguments. Please see help file for usage');
elseif nargin > 10
    error('Too many input arguments.  Please see help file for usage');
else
    PARAM.filename=varargin{1};
    PARAM.loadevents.state=varargin{2};
    PARAM.loadevents.type=varargin{3};
    PARAM.loadevents.dig_ch1=varargin{4};
    PARAM.loadevents.dig_ch1_label=varargin{5};
    PARAM.loadevents.dig_ch2=varargin{6};
    PARAM.loadevents.dig_ch2_label=varargin{7};
    PARAM.chan_adjust_status=varargin{8};
    PARAM.chan_adjust=varargin{9};
    PARAM.chans=varargin{10};
end
Trigs1=[];
Trigs2=[];



% ---------------- Opening File------------------
trcfile=PARAM.filename;
fid=fopen(trcfile,'r','n','windows-1252');
if fid==-1
    error('Can''t open *.trc file')
end


%------------------reading patient & recording info----------
fseek(fid,64,-1);
surname=char(fread(fid,22,'char'))';
name=char(fread(fid,20,'char'))';

fseek(fid,128,-1);
day=fread(fid,1,'char');
if length(num2str(day))<2
    day=['0' num2str(day)];
else
    day=num2str(day);
end
month=fread(fid,1,'char');
switch month
case 1 
    month='JAN';
case 2 
    month='FEB';
case 3 
    month='MAR';
case 4 
    month='APR';
case 5 
    month='MAY';
case 6 
    month='JUN';
case 7 
    month='JUL';
case 8 
    month='AUG';
case 9 
    month='SEP';
case 10 
    month='OCT';
case 11 
    month='NOV';
case 12 
    month='DEC';
end
year=num2str(fread(fid,1,'char')+1900);

%------------------ Reading Header Info ---------

fseek(fid,175,-1);
Header_Type=fread(fid,1,'char');
if Header_Type ~= 4
    error('*.trc file is not Micromed System98 Header type 4')
end

fseek(fid,138,-1);
Data_Start_Offset=fread(fid,1,'uint32');
Num_Chan=fread(fid,1,'uint16');
Multiplexer=fread(fid,1,'uint16');
Rate_Min=fread(fid,1,'uint16');
Bytes=fread(fid,1,'uint16');
fseek(fid,176+8,-1);
Code_Area=fread(fid,1,'uint32');
Code_Area_Length=fread(fid,1,'uint32');
fseek(fid,192+8,-1);
Electrode_Area=fread(fid,1,'uint32');
Electrode_Area_Length=fread(fid,1,'uint32');

fseek(fid,400+8,-1);
Trigger_Area=fread(fid,1,'uint32');
Tigger_Area_Length=fread(fid,1,'uint32');


%----------------- Allocate Memory and Determine Data Type----------

fseek(fid,Data_Start_Offset,-1);

fprintf('Allocating memory...\n');
switch Bytes
case 1   
    bstring='uint8';
case 2
    bstring='uint16';
case 4
    bstring='uint32';
end
       
tracetmp=fread(fid,bstring,(Num_Chan-1)*Bytes)';
traceL=length(tracetmp);
clear tracetmp;

if isempty(PARAM.chans) 
    chans=1:Num_Chan;
else
    chans=eval([ '[' PARAM.chans ']' ]);
end

chansL=length(chans);
tracedata=zeros(chansL,traceL);
m=traceL;

%------------------ Reading Code Info -------------
fseek(fid,Code_Area,-1);
code=fread(fid,Num_Chan,'uint16');


for c=1:Num_Chan
    electrode(c).chan_record=code(c);
    fseek(fid,Electrode_Area+code(c)*128,-1);
    fseek(fid,2,0);
    if c <10 
        electrode(c).positive_input=...
            [num2str(c),' -',char(fread(fid,6,'char'))'];
        electrode(c).negative_input=...
            [num2str(c),' -',char(fread(fid,6,'char'))'];
    else
        electrode(c).positive_input=...
            [num2str(c),'-',char(fread(fid,6,'char'))'];
        electrode(c).negative_input=...
            [num2str(c),'-',char(fread(fid,6,'char'))'];
    end
    electrode(c).logical_min=fread(fid,1,'int32');
    electrode(c).logical_max=fread(fid,1,'int32');
    electrode(c).logical_ground=fread(fid,1,'int32');
    electrode(c).physical_min=fread(fid,1,'int32');
    electrode(c).physical_max=fread(fid,1,'int32');
    
    electrode(c).measurement_unit=fread(fid,1,'int16');
    switch electrode(c).measurement_unit
    case -1
        electrode(c).measurement_unit=1e-9; 
    case 0
        electrode(c).measurement_unit=1e-6;
    case 1
        electrode(c).measurement_unit=1e-3;
    case 2
        electrode(c).measurement_unit=1;
    case 100
        electrode(c).measurement_unit='percent';
    case 101
        electrode(c).measurement_unit='bpm';
    case 102
        electrode(c).measurement_unit='Adim';
    otherwise
        warning('Unknown measurement unit. uV assumed.');
        electrode(c).measurement_unit=10e-6;
    end
    fseek(fid,8,0);
    electrode(c).rate_coef=fread(fid,1,'uint16'); 
end


%---------------- Read & Prep Trigger Area Data ----------
fseek(fid,Trigger_Area,-1);
for l=1:Tigger_Area_Length/6
    trigger(1,l)=fread(fid,1,'uint32');
    trigger(2,l)=fread(fid,1,'uint16');
end


first_trigger=trigger(1,1);
tl=length(trigger);
NoTrig=0;
for tr=1:tl
    if ((trigger(1,tr) <= m) & (trigger(1,tr) >= first_trigger))
        NoTrig=NoTrig+1;
    end
end

if NoTrig > 0
   	trigger=trigger(:,1:NoTrig);
else
	trigger=[];
	first_trigger=[];
end

%---------------Reading Other Event Data   -------------

switch  PARAM.loadevents.state
case 'no';
case 'yes';
    fprintf('Extracting events...\n');
    switch PARAM.loadevents.type
    case 'marker'
        if ~isempty(trigger)
            [triggerR,triggerC]=size(trigger)
            for E=1:triggerC
                TRC.event(end+1).type=num2str(trigger(2,E));
                TRC.event(end).latency=trigger(1,E)+1;
            end
        else
            warndlg('No marker triggers to import on ''MRK'' channel',...
                'Import .TRC Warning!');
        end
    case 'eegchan'
        %------find 1st trigger channel-------
        dig_ch=[];
        if str2num(PARAM.loadevents.dig_ch1)
            dig_ch(1)=str2num(PARAM.loadevents.dig_ch1);
        else
            for C=1:Num_Chan
                if ~isempty(findstr(lower(PARAM.loadevents.dig_ch1),...
                        lower(electrode(C).positive_input)))
                    dig_ch(1)=C;
                    break;
                end
            end
        end
        
        %------check for 2nd trigger channel-------
        if ~isempty(PARAM.loadevents.dig_ch2)
            if str2num(PARAM.loadevents.dig_ch2)
                dig_ch(2)=str2num(PARAM.loadevents.dig_ch2);
            else
                for C=1:Num_Chan
                    if ~isempty(findstr(lower(PARAM.loadevents.dig_ch2),...
                            lower(electrode(C).positive_input)))
                        dig_ch(2)=C;
                        break;
                    end
                end
            end
        end
        
                
        %--------read trigs-------------------------
  
        fseek(fid,Data_Start_Offset+(dig_ch(1)-1)*Bytes,-1);
        tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
        if ischar(electrode(dig_ch(1)).measurement_unit)==0
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min)*...
                electrode(dig_ch(1)).measurement_unit;
        else
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min);
        end
        trigval=min(tracedata(1,:));
        Trigs1=find(tracedata(1,:)==trigval);
        

        
        if length(dig_ch)>1
            fseek(fid,Data_Start_Offset+(dig_ch(2)-1)*Bytes,-1);
            tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
            if ischar(electrode(dig_ch(2)).measurement_unit)==0
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min)*...
                    electrode(dig_ch(2)).measurement_unit;
            else
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min);
            end
            Trigs2=find(tracedata(1,:)==trigval);
        end
        
      
        %--------write trigs----------------------------------
        for E=1:length(Trigs1)
            TRC.event(end+1).type=PARAM.loadevents.dig_ch1_label;
            TRC.event(end).latency=Trigs1(E);
        end
            
        if ~isempty(Trigs2)
            for E=1:length(Trigs2)
                TRC.event(end+1).type=PARAM.loadevents.dig_ch2_label;
                TRC.event(end).latency=Trigs2(E);
            end
          
        end
        
    case 'both'
        
        %-----------Marker trigs----------------------------
        if ~isempty(trigger)
            [triggerR,triggerC]=size(trigger);
            for E=1:triggerC
                TRC.event(end+1).type=num2str(trigger(2,E));
                TRC.event(end).latency=trigger(1,E)+1;
            end
        else
            warndlg('No marker triggers to import on ''MRK'' channel',...
                'Import .TRC Warning!');
        end
        
         %------find 1st trigger channel-------
        dig_ch=[];
        if str2num(PARAM.loadevents.dig_ch1)
            dig_ch(1)=str2num(PARAM.loadevents.dig_ch1);
        else
            for C=1:Num_Chan
                if ~isempty(findstr(lower(PARAM.loadevents.dig_ch1),...
                        lower(electrode(C).positive_input)))
                    dig_ch(1)=C;
                    break;
                end
            end
        end
        
        %------check for 2nd trigger channel-------
        if ~isempty(PARAM.loadevents.dig_ch2)
            if str2num(PARAM.loadevents.dig_ch2)
                dig_ch(2)=str2num(PARAM.loadevents.dig_ch2);
            else
                for C=1:Num_Chan
                    if ~isempty(findstr(lower(PARAM.loadevents.dig_ch2),...
                            lower(electrode(C).positive_input)))
                        dig_ch(2)=C;
                        break;
                    end
                end
            end
        end
        
                
        %--------read trigs-------------------------
       
        fseek(fid,Data_Start_Offset+(dig_ch(1)-1)*Bytes,-1);
        tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
        if ischar(electrode(dig_ch(1)).measurement_unit)==0
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min)*...
                electrode(dig_ch(1)).measurement_unit;
        else
            tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(1)).logical_ground)/...
                (electrode(dig_ch(1)).logical_max-...
                electrode(dig_ch(1)).logical_min+1))*...
                (electrode(dig_ch(1)).physical_max-...
                electrode(dig_ch(1)).physical_min);
        end
        trigval=min(tracedata(1,:));
        Trigs1=find(tracedata(1,:)==trigval);
        

        
        if length(dig_ch)>1
            fseek(fid,Data_Start_Offset+(dig_ch(2)-1)*Bytes,-1);
            tracedata(1,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
            if ischar(electrode(dig_ch(2)).measurement_unit)==0
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min)*...
                    electrode(dig_ch(2)).measurement_unit;
            else
                tracedata(1,:)=-((tracedata(1,:)-electrode(dig_ch(2)).logical_ground)/...
                    (electrode(dig_ch(2)).logical_max-...
                    electrode(dig_ch(2)).logical_min+1))*...
                    (electrode(dig_ch(2)).physical_max-...
                    electrode(dig_ch(2)).physical_min);
            end
            Trigs2=find(tracedata(1,:)==trigval);
        end
        
        %------------write trigs----------------------------------
        
        for E=1:length(Trigs1)
            TRC.event(end+1).type=PARAM.loadevents.dig_ch1_label;
            TRC.event(end).latency=Trigs1(E);
        end
            
        if ~isempty(Trigs2)
            for E=1:length(Trigs2)
                TRC.event(end+1).type=PARAM.loadevents.dig_ch2_label;
                TRC.event(end).latency=Trigs2(E);
            end
          
        end
        
    end   
end


%------------------Reading Data-------------------

fprintf('Reading data...\n');
for c=1:chansL
    fseek(fid,Data_Start_Offset+(chans(c)-1)*Bytes,-1);
    tracedata(c,:)=fread(fid,bstring,(Num_Chan-1)*Bytes)';
    if ischar(electrode(chans(c)).measurement_unit)==0
        tracedata(c,:)=-((tracedata(c,:)-electrode(chans(c)).logical_ground)/...
            (electrode(chans(c)).logical_max-...
            electrode(chans(c)).logical_min+1))*...
            (electrode(chans(c)).physical_max-...
            electrode(chans(c)).physical_min)*...
            electrode(chans(c)).measurement_unit;
    else
        tracedata(c,:)=-((tracedata(c,:)-electrode(chans(c)).logical_ground)/...
            (electrode(chans(c)).logical_max-...
            electrode(chans(c)).logical_min+1))*...
            (electrode(chans(c)).physical_max-...
            electrode(chans(c)).physical_min);
    end
    
    if ~isempty(Trigs1)
        if chans(c)==dig_ch(1)
            tracedata(c,Trigs1)=...
                (tracedata(c,(Trigs1+1))+tracedata(c,(Trigs1-1)))/2;
        end
    end

    if ~isempty(Trigs2)
        if chans(c)==dig_ch(2)
            tracedata(c,Trigs2)=...
                (tracedata(c,(Trigs2+1))+tracedata(c,(Trigs2-1)))/2;
        end
    end   
    
end


% -----------Reading Fs-------------------------

mean_fs=mean(cat(1,electrode.rate_coef));
switch mean_fs
case 1
    fs=1*Rate_Min;
case 2
    fs=2*Rate_Min;
case 3
    fs=3*Rate_Min;
case 4
    fs=4*Rate_Min;
case 5
    fs=5*Rate_Min;
otherwise
    warning('Unsupported Sampling Frequency');
end


%----------Prep output-------------------------

fprintf('Preparing output...\n');
if PARAM.chan_adjust_status==1
	if length(tracedata)>(fs*60)
        avgvar=mean(var(tracedata(1:Num_Chan-3,1:fs*60)'));
	else
        avgvar=mean(var(tracedata(1:Num_Chan-3,:)'));
	end
    
    ch_adj_t=['[' PARAM.chan_adjust ']'];
    ch_adj=eval(ch_adj_t);
    
    for ch=1:length(ch_adj)
        tracedata(ch_adj(ch),:)=tracedata(ch_adj(ch),:)*avgvar/var(tracedata(ch_adj(ch),:));
    end
end


TRC.data=-tracedata*1e6; % scale to uV and change polarity for EEGLAB
sp=findstr(surname,'  ');
if sp >=1
    TRC.setname=[surname(1:(sp-1)) ', ' name(1) '. ' year month day ' .TRC File'];
else
    TRC.setname=[surname ', ' name(1) '. ' year month day ' .TRC File'];
end
TRC.filename=trcfile;
TRC.filepath='';
TRC.pnts=length(tracedata);
TRC.nbchan=chansL;
TRC.trials=1;
TRC.srate=fs;
TRC.xmin=0;
TRC.xmax=(TRC.pnts-1)/fs;
TRC.ref='common';
return;



