% read_edf() - read eeg data in EDF+ format.
%
% Orignal Code-file:
% Jeng-Ren Duann, CNL/Salk Inst., 2001-12-21
%
% Modifications:
% 03-21-02 editing hdr, add help -ad 
% 05-03-12 FHatz Neurology Basel (Support for edf+)
% 05-07-12 A. Barborica Bucharest University - fixed error when
% variable sample rate was used for different channels
% 05-08-12 C. Donos - implemented reading of a specified interval from the
% entire recording
% 06-11-09 O. David - added as is in ImaGIN toolbox

% Usage: 
%    >> [data,header] = read_edf(filename,tstart,tstop,voltage_sw);
%
% Input:
%    filename - file name of the eeg data
%    tstart   - start time, in seconds, of the segment to be loaded
%    tstop    - end time, in seconds, of the segment to be loaded
%    voltage_sw - switch for returning physical values instead of digital conversion values
%
% 
% Output:
%    data   - eeg data in ADC samples (channel, timepoint)
%    header - structured information about the read eeg data
%      header.events - events (structure: .POS .DUR .TYP)
%      header.numtimeframes - length of EEG data
%      header.samplingrate = samplingrate
%      header.numchannels - number of channels
%      header.numauxchannels - number of non EEG channels (only ECG channel is recognized)
%      header.channels - channel labels
%      header.year - timestamp recording
%      header.month - timestamp recording
%      header.day - timestamp recording
%      header.hour - timestamp recording
%      header.minute - timestamp recording
%      header.second - timestamp recording
%      header.ID - EEG number
%      header.technician - responsible investigator or technician
%      header.equipment - used equipment
%      header.subject.ID - local patient identification
%      header.subject.sex - M or F
%      header.subject.name - patients name
%      header.subject.year - birthdate
%      header.subject.month - birthdate
%      header.subject.day - birthdate
%      header.hdr - original header

function [data,header] = ImaGIN_read_edf(filename,tstart,tstop,voltage_sw)


if nargin < 1
    help readedf;
    return;
end;
    
fp = fopen(filename,'r','ieee-le');
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

if (nargin < 4)
    voltage_sw = true;
end
if (nargin < 3)
    recstop = hdr.records;
else
    recstop = min(hdr.records,round(tstop/hdr.duration));
end
if (nargin < 2)
    recstart = 0;
else
    recstart = min(recstop,round(tstart/hdr.duration));
end

offset=hdr.length + recstart*hdr.duration*sum(hdr.numbersperrecord)*2; %% 16 bit

fseek(fp,offset,-1);

% fseek(fp,hdr.length,-1);
% data = fread(fp,'int16');
% data = fread(fp,(recstop-recstart)*hdr.duration*sum(hdr.numbersperrecord),'int16');
data = fread(fp, (recstop-recstart)*sum(hdr.numbersperrecord), 'int16');   % FIX FT 14-Mar-2017
%fclose(fp);

header.hdr = hdr;
%header.samplingrate = hdr.numbersperrecord(1) / hdr.duration;
header.samplingrate = hdr.numbersperrecord / hdr.duration;  % AB: Previous formula was incorrect, as different channels may have different sampling rates.
header.numchannels = hdr.channels;
header.numauxchannels = 0;
header.channels = hdr.channelname;
tmp = textscan(hdr.intro(89:168),'%s');
tmp = tmp{1,1};
[header.year header.month header.day] = datevec(tmp(2,1));
header.hour = str2num(hdr.intro(177:178));
header.minute = str2num(hdr.intro(180:181));
header.second = str2num(hdr.intro(183:184));
try
    header.startdatetime = datestr(datenum(hdr.intro(169:184),'dd.mm.yyHH.MM.SS'));
end
header.ID = tmp{3,1};
header.technician = tmp{4,1};
header.equipment = tmp{5,1};
tmp = textscan(hdr.intro(9:88),'%s');
tmp = tmp{1,1};
if (size(tmp,1)>=1)
    header.subject.ID = tmp{1,1};
else
    header.subject.ID = 'X';
end
if (size(tmp,1)>=2)
    header.subject.sex = tmp{2,1};
else
    header.subject.sex = 'X';
end
if (size(tmp,1)>=4)
    header.subject.name = tmp{4,1};
else
    header.subject.name = 'X';
end
if (size(tmp,1)>=3) && (~strcmp(tmp(3,1),'X'))
    [header.subject.year header.subject.month header.subject.day] = datevec(tmp(3,1));
else
    [header.subject.year header.subject.month header.subject.day] = datevec('Jan 01 1900');
end
%clearvars tmp
clear tmp

if strcmp(hdr.channelname(end,1:15),'EDF Annotations') || strcmp(hdr.channelname(end,1:6),'Status')
    ix=find(hdr.numbersperrecord==hdr.numbersperrecord(1));   % AB: Find all channels that have the same sample rate as the first channel
    hix=1:hdr.numbersperrecord(1);  %OD
    for i0=2:length(ix)
        hix=[hix sum(hdr.numbersperrecord(1:ix(i0)-1))+[1:hdr.numbersperrecord(1)]];
    end
    n = length(hdr.numbersperrecord) - ix(end);               % AB: Remove last n channels %Corrected by OD
    n = length(hdr.numbersperrecord) - length(ix);               % OD
    header.numchannels = header.numchannels - n;
%     header.channels = header.channels(1:end-n,:);
    header.channels = header.channels(ix,:);    %OD
    try
        data = reshape(data,sum(hdr.numbersperrecord),(recstop-recstart));
    catch
        warning('Size of data does not match the header, ugly fix used by OD')
        tmp=floor(length(data)/sum(hdr.numbersperrecord));
        recstop=tmp+recstart;
        data = reshape(data(1:sum(hdr.numbersperrecord)*tmp),sum(hdr.numbersperrecord),(recstop-recstart));
    end
    eventstmp = data((end - hdr.numbersperrecord(end) + 1):end,:);
%     data = data(1:end-sum(hdr.numbersperrecord((end-n+1):end)),:);
    data = data(hix,:); %OD
    data = reshape(data,hdr.numbersperrecord(1),header.numchannels,(recstop-recstart));
    data = permute(data,[1 3 2]);
    data = reshape(data,hdr.numbersperrecord(1)*(recstop-recstart),header.numchannels);
    % Read annotations from beginning file
    fseek(fp,hdr.length,-1);
    i = 0;
    recdata = fread(fp,2*hdr.duration*sum(hdr.numbersperrecord),'int8');
    evttmp = recdata((end - 2*hdr.numbersperrecord(end) + 1):end);
    while (evttmp(15) > 0)
        i = i + 1; 
        eventsall(i,:) = evttmp;
        % Read next record
        recdata = fread(fp,2*hdr.duration*sum(hdr.numbersperrecord),'int8');
        evttmp = recdata((end - 2*hdr.numbersperrecord(end) + 1):end);
    end
    % Read additional events that may be saved in the specified segment
    for j = (max(recstart,i)+1 - recstart):(recstop-recstart)
        evttmp = typecast(int16(eventstmp(:,j)),'uint8')';
        if (evttmp(15) > 0)
            i = i + 1; 
            eventsall(i,:) = evttmp;
        end
    end
    header.events.TYP = [];
    header.events.POS = [];
    header.events.TIME = [];
    header.events.DUR = [];
    if exist('eventsall','var')
        eventsmod = eventsall;
        eventsmod(find(eventsall == 32)) = 95;
        eventsmod(find(eventsall == 20)) = 32;
        eventsmod(find(eventsall == 21)) = 32;
        eventsmod(find(eventsall == 43)) = 32;
        eventsmod(find(eventsall == 0)) = 32;
        for i = 1:size(eventsall,1)
            eventspos = find(eventsall(i,:) == 43);
            for j = 2:length(eventspos)
                try %OD
                    tmp = textscan(native2unicode(eventsmod(i,eventspos(j):end)),'%s');
                    tmp = tmp{1,1};
                    evtpos = str2double(tmp(1,1));
                    if (evtpos>=recstart) && (evtpos<=recstop)
                        header.events.POS  = [header.events.POS  (str2double(tmp(1,1))*max(header.samplingrate))];
                        header.events.TIME = [header.events.TIME (str2double(tmp(1,1)))];
                        if (size(tmp,1) > 1)
%                             if str2double(tmp(2,1)) > 0
                            if str2double(tmp(2,1)) >= 0        %OD
                                header.events.DUR = [header.events.DUR (str2double(tmp(2,1))*max(header.samplingrate))];
                                if size(tmp,1) > 2
                                    header.events.TYP = [header.events.TYP tmp(3,1)];
                                else
                                    header.events.TYP = [header.events.TYP tmp(2,1)];
                                end
                            else
                                header.events.DUR = [header.events.DUR 1];
                                header.events.TYP = [header.events.TYP tmp(2,1)];
                            end
                        else
                            header.events.DUR = [header.events.DUR 1];
                            header.events.TYP = [header.events.TYP ' '];
                        end
                    end
                end
            end
        end
    end
    
else
    data = reshape(data,hdr.numbersperrecord(1),hdr.channels,(recstop-recstart));
    temp = [];
    for i=1:(recstop-recstart),
        temp = [temp data(:,:,i)'];
    end
    data = temp;
end
fclose(fp);

if strcmp(header.channels(end,1:3),'ECG')
    header.numauxchannels = 1;
    header.numdatachannels = header.numchannels -1;
else
    header.numdatachannels = header.numchannels;
end
header.numtimeframes = size(data,2);

if (voltage_sw)
    for i=1:size(data,2)
           data(:,i) = data(:,i) * (hdr.physmax(i)-hdr.physmin(i))/(hdr.digimax(i)-hdr.digimin(i));
    end
end