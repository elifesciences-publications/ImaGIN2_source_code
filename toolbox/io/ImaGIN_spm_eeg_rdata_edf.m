function D = ImaGIN_spm_eeg_rdata_edf(S)
% converts EEG data from edf - to SPM-format with monopolar montage for SEEG
% FORMAT D = ImaGIN_spm_eeg_rdata_edf(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of Deltamed-file
% Fchannels   - filename of channel template
% reference   - name of reference channel(s)
% channel     - index of channel(s) to read
%_______________________________________________________________________
% 
% ImaGIN_spm_eeg_rdata_deltamedbin reads a *.bin Deltamed file and its associated header (.txt), stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
%_______________________________________________________________________

% Olivier David


D = [];


try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.edf$', 'Select bin file');
end
cd(spm_str_manip(Fdata,'h'))
FdataName=spm_str_manip(Fdata,'t');
FdataName=FdataName(1:end-4);

try
    FileOut=S.FileOut;
catch
    FileOut=FdataName;
end

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
    S.SizeMax;
catch
    S.SizeMax = str2num(spm_input('Maximum samples per file ', '+1', 's',1e12));
end

try
    S.SEEG;
catch
    S.SEEG = spm_input('SEEG? ','+1','Yes|No');
end
switch S.SEEG
    case 'Yes'
        S.SEEG=1;
    case 'No'
        S.SEEG=0;
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



%Estimates the number of files
fp = fopen(Fdata,'r','ieee-le');
if fp == -1,
  error('File not found ...!');
  return;
end
hdr.intro = setstr(fread(fp,256,'uchar')');
hdr.startdate = hdr.intro(169:176);
hdr.starttime = hdr.intro(177:184);
try
    hdr.startdatetime = datestr(datenum(hdr.intro(169:184),'dd.mm.yyHH.MM.SS'));
end
hdr.length = str2num(hdr.intro(185:192));
hdr.records = str2num(hdr.intro(237:244));
hdr.duration = str2num(hdr.intro(245:252));
hdr.channels = str2num(hdr.intro(253:256));
hdr.channelname = setstr(fread(fp,[16,hdr.channels],'char')');
hdr.transducer = setstr(fread(fp,[80,hdr.channels],'char')');
hdr.physdime = setstr(fread(fp,[8,hdr.channels],'char')');
hdr.physmin = str2num(setstr(fread(fp,[8,hdr.channels],'char')'));
hdr.physmax = str2num(setstr(fread(fp,[8,hdr.channels],'char')'));
hdr.digimin = str2num(setstr(fread(fp,[8,hdr.channels],'char')'));
hdr.digimax = str2num(setstr(fread(fp,[8,hdr.channels],'char')'));
hdr.prefilt = setstr(fread(fp,[80,hdr.channels],'char')');
hdr.numbersperrecord = str2num(setstr(fread(fp,[8,hdr.channels],'char')'));
fclose(fp)
ndata=sum(hdr.numbersperrecord)*hdr.records;
if ndata<S.SizeMax
    NSegments=1;
    tstart=0;
    tend=inf;
else
    NSegments=ceil(ndata/S.SizeMax);
    L=round(hdr.records/NSegments);
    tstart=0:L:hdr.records-1;
    tstart=tstart(1:NSegments);
    tstart=round(tstart-0.02*L);
    tstart(find(tstart<0))=0;
    tend=L:L:hdr.records+L;
    tend=tend(1:NSegments);
    tend=round(tend+0.02*L);
%     tend(find(tend>hdr.hdr.records))=hdr.hdr.records;
end

% [Data hdr] = ImaGIN_read_edf(Fdata,1,2);
% ndata=sum(hdr.numbersperrecord)*hdr.records;
% if ndata<S.SizeMax
%     NSegments=1;
%     tstart=0;
%     tend=inf;
% else
%     NSegments=ceil(ndata/S.SizeMax);
%     L=round(hdr.hdr.records/NSegments);
%     tstart=1:L:hdr.hdr.records;
%     tstart=tstart(1:NSegments);
%     tstart=round(tstart-0.02*L);
%     tstart(find(tstart<0))=0;
%     tend=L:L:hdr.hdr.records+L;
%     tend=tend(1:NSegments);
%     tend=round(tend+0.02*L);
% %     tend(find(tend>hdr.hdr.records))=hdr.hdr.records;
% end

Schannel=S.channel;
for i00=1:NSegments
    
    %Read edf
    [Data hdr] = ImaGIN_read_edf(Fdata,tstart(i00),tend(i00));
    Data=transpose(Data);
    
    if isempty(Schannel)
        %remove channels that do not have the same sampling rate as the first one
        Remove=[];
        for i1=2:hdr.numchannels
            if hdr.samplingrate(i1)~=hdr.samplingrate(1)
                Remove(end+1)=i1;
            end
        end
        Keep=setdiff(1:hdr.numchannels,Remove);
    else
        Keep=Schannel;
    end
    Data=Data(Keep,:);
    hdr.numchannels=length(Keep);
    hdr.numtimeframes=length(Keep);
    hdr.numdatachannels=length(Keep);
    hdr.channels=hdr.channels(Keep,:);
    S.channel=1:length(Keep);
    
    
    
    if S.coarse>1
        for i1=1:size(Data,1)
            Data(i1,:)= ImaGIN_lowpassFilter(Data(i1,:),hdr.samplingrate(i1),hdr.samplingrate(i1)/(2*S.coarse)-2);    %antialiasing,
        end
        Data=Data(:,1:S.coarse:end);
    end
    
    
    
    %save as SPM file
    D = [];
    D.Fsample=hdr.samplingrate(1)/S.coarse;
    Nsamples = size(Data,2);
    Nchannels = size(Data,1);
    D.channels = repmat(struct('bad', 0), Nchannels,1);
    
    D.timeOnset=tstart(i00);
    
    %events
    evt = [];
    for i0=1:length(hdr.events.TYP)
        evt(i0).type  = hdr.events.TYP{i0};
        evt(i0).time = hdr.events.TIME(i0);
        evt(i0).value= hdr.events.TYP{i0};
    end
    
    if isfield(D,'time')
        time=D.time;
    else
        time=[0:1/D.Fsample:(Nsamples-1)/D.Fsample];
        time=time+D.timeOnset;
    end
    
    
    % Read the electrode names
    Name={};
    for i1=1:size(Data,1)
        
        Name{i1,1}=lower(deblank(hdr.channels(i1,:)));
        
        %remove zeros
        a=Name{i1,1};
        v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
        if length(v_num)>1
            if v_num(2)==v_num(1)+1 && strcmp(a(v_num(1)),'0')
                Name{i1,1}=a(setdiff(1:length(a),v_num(1)));
            end
        end
        
    end
    
    %Assign types to electrodes
    SEEGMNI=0;
    isecg=[];
    KeepMNI=[];
    for i1=1:length(S.channel)
        D.channels(i1).label=Name(S.channel(i1));
        switch D.channels(i1).label{1}
            case{'ecg1','ecg2','ekg1','ekg2'}
                D.channels(i1).type='ECG';
                isecg(end+1)=i1;
            otherwise
                D.channels(i1).type='EEG';
        end
        
        %clean label (for MNI)
        [tmp,tmp2]=strtok(Name{i1},' ');
        if i1==1&&S.SEEG
            if strcmp(tmp,'seeg')
                SEEGMNI=1;
            end
        end
        if SEEGMNI
            switch tmp
                case 'seeg'
                    KeepMNI(end+1)=i1;
                    D.channels(i1).type='EEG';
                    [tmp1,tmp2]=strtok(strtrim(tmp2),'-');
                    D.channels(i1).label=tmp1;
                case 'eeg'
                    D.channels(i1).type='EEG';
                    D.channels(i1).label=Name{i1};
                case {'ekg','ecg'}
                    KeepMNI(end+1)=i1;
                    D.channels(i1).type='ECG';
                    [tmp1,tmp2]=strtok(strtrim(tmp2),'-');
                    D.channels(i1).label=tmp1;
                    isecg=[isecg length(KeepMNI)];
            end
        else
            switch tmp
                case 'eeg'
                    D.channels(i1).type='EEG';
                    [tmp1,tmp2]=strtok(strtrim(tmp2),'-');
                    switch tmp2
                        case '-ref'
                            D.channels(i1).label=tmp1;
                        otherwise
                            D.channels(i1).label=Name{i1};
                    end
                case {'ekg','ecg'}
                    D.channels(i1).type='ECG';
                    [tmp1,tmp2]=strtok(strtrim(tmp2),'-');
                    switch tmp2
                        case '-ref'
                            D.channels(i1).label=tmp1;
                        otherwise
                            D.channels(i1).label=Name{i1};
                    end
                    isecg=[isecg i1];
                otherwise
                            D.channels(i1).label=Name{i1};
                    
            end
        end
    end
    %remove _
    for i1=1:length(S.channel)
    	[tmp1,tmp2]=strtok(D.channels(i1).label,'_');
        if length(tmp2)>1
            D.channels(i1).label=[tmp1 tmp2(2:end)];
        end
    end
        

    if SEEGMNI
        Data=Data(KeepMNI,:);
        S.channel=1:length(KeepMNI);
        Nchannels=length(KeepMNI);
        D.channels=D.channels(KeepMNI);
        if ~isempty(isecg)
            if length(intersect(KeepMNI,isecg))>0
                isecg=length(find(KeepMNI)<=isecg);
            end
        end
    end
    

    
    
    %reorder to put ECG channels at the end
    reorder=[setdiff(1:length(S.channel),isecg) isecg];
    D.channels=D.channels(reorder);
    Data=Data(reorder,:);
    
    
    %Remove non SEEG channels
    if S.SEEG
        Keep=[];
        for i1=1:length(D.channels)
            s_test=1;
            a_string=D.channels(i1).label;
            a=a_string;
            v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
            
            
            %assumes all SEEG channels are the first in the file, their last
            %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter
            if (isempty(v_num))
                s_test=0;
            elseif (length(v_num)>2)
                s_test=0;
            elseif (max(diff(v_num))>1)
                s_test=0;
            else
                s_test=((max(v_num)-min(v_num))==(length(v_num)-1));
                s_test=s_test*(max(v_num)==length(a_string));
                s_test=s_test*(isletter(a_string(1)));
                a_rt=a_string(1:(min(v_num)-1));
                if strcmp(lower(a_rt),'myo')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'oc')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'dd')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'dc')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'dg')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'el')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'eog')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'myog')
                    s_test=0;
                end
                if strcmp(lower(a_rt),'myod')
                    s_test=0;
                end
                try
                    if strcmp(lower(a_rt(1:3)),'emg')
                        s_test=0;
                    end
                end
                try
                    if strcmp(lower(a_rt(1:5)),'photic')
                        s_test=0;
                    end
                end
                
            end
            if s_test
                Keep=[Keep i1];
            end
               
        end
        Data=Data(Keep,:);
        S.channel=1:length(Keep);
        Nchannels=length(Keep);
        D.channels=D.channels(Keep);
        if ~isempty(isecg)
            if length(intersect(Keep,isecg))>0
                isecg=length(find(Keep)<=isecg);
            end
        end

    end

    
    
    D.Nchannels = Nchannels;
    D.Nsamples = Nsamples;
    nsampl=Nsamples;
    nchan=Nchannels;
    
    %Assume continuous data
    D.trials.label = 'Undefined';
    D.trials.events = [];
    D.trials.onset = 1/D.Fsample;
    
    
    
    %Continuous
    if NSegments==1
        D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
    else
        D.data = file_array([FileOut '_' num2str(i00) '.dat'], [nchan nsampl], 'float32-le');
    end
    % physically initialise file
    D.data(end,end) = 0;
    offset = 1;
    nblocksamples = size(Data,2);
    D.data(:, offset:(offset+nblocksamples-1)) = full(Data);
    
    
    %Electrodes
    for i1=1:length(D.channels)
        if iscell(D.channels(i1).label)
            D.sensors.eeg.label(i1)=D.channels(i1).label;
        else
            D.sensors.eeg.label{1,i1}=D.channels(i1).label;
        end
    end
    D.sensors.eeg.unit='mm';
    D.sensors.eeg.pnt=NaN*zeros(length(D.channels),3);
    D.sensors.eeg.type='eeg';
    D.sensors.chantype='eeg';
    
    D.Atlas=S.Atlas;
    
    %--------- Create meeg object
    if NSegments==1
        [DirOut,NameOut,~]=fileparts(FileOut);
    else
        [DirOut,NameOut,~]=fileparts([FileOut '_' num2str(i00)]);
    end
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
    
    
%     %clean events
%     ev=events(D);
%     Select=[];
%     for i1=1:length(ev)
%         if ~isempty(ev(i1).value)&&~strcmp(ev(i1).value(1),'+')
%             Select=[Select i1];
%         end
%     end
%     for i1=1:length(Select)
%         ev(Select(i1)).type=ev(Select(i1)).value;
%     end
%     if ~isempty(Select)
%         D=events(D,1,ev(Select));
%     end
%     save(D);
        
    
end
spm('Pointer','Arrow');
    
