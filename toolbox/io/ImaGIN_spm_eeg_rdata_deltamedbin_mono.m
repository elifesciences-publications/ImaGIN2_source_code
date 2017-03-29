function D = ImaGIN_spm_eeg_rdata_deltamedbin_mono(S)
% Converts EEG data from binary Deltamed (.bin and .txt)- to SPM-format
% 
% USAGE:   D = ImaGIN_spm_eeg_rdata_deltamedbin(S)

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

D = [];


try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.bin$', 'Select bin file');
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

%Read header
Fdatatxt=[Fdata(1:end-4) '.txt'];
header=readheaderbin(Fdatatxt);

% Intracranial human epi
try
    S.SEEG;
catch
    S.SEEG = spm_input('SEEG? ','+1','Yes|No');
end
switch S.SEEG
    case{'Yes'}
        v_rt=[];
        s_counter=1;
        s_other=1;
        v_neworder=[];
        v_name=[];
        m_bipole_ECG=[];
        EEGFlag=1;
        for s_c=1:length(header.ChannelNames)
            a_string=header.ChannelNames{s_c};
            a_string=deblank(a_string);
            
            %OD for Lyon
            a_string=lower(a_string(~isspace(a_string)));
            %OD for Brousse
            a_string=a_string(~(a_string=='.'));
            ss_channel.values(s_c).value=a_string;
            
            a=a_string;
            v_name=strvcat(v_name,a_string);
            v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
            
            
            %assumes all SEEG channels are the first in the file, their last
            %characters are 1 or 2 numbers (no more), no other numbers in the file name, first character is a letter
            if (isempty(v_num))
                s_test=0;
                EEGFlag=0;
                a_rt=a_string;
                if strcmp(lower(a_rt),'ecg')
                    s_test=1;
                    EEGFlag=1;
                    m_bipole_ECG=[s_c s_c];
                end
            elseif (length(v_num)>2)
                s_test=0;
                EEGFlag=0;
            elseif (max(diff(v_num))>1)
                s_test=0;
                EEGFlag=0;
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
%                     EEGFlag=0;
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
        m_bipole=[m_bipole;m_bipole_ECG];
        
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
        
        f=fopen(fullfile(spm_str_manip(Fdatatxt,'h'),'bipole.txt'),'w');
        for s_b=1:size(m_bipole,1)
%             fprintf(f,'%s\t%s\n',ss_channel.values(v_neworder_rev(m_bipole(s_b,1))).value,ss_channel.values(v_neworder_rev(m_bipole(s_b,2))).value);
            fprintf(f,'%i\t%i\n',v_neworder_rev(m_bipole(s_b,1)),v_neworder_rev(m_bipole(s_b,2)));
        end;
        fclose(f);
        S.bipole=load(fullfile(spm_str_manip(Fdatatxt,'h'),'bipole.txt'));
        S.Bipolar='Yes';
        
        
        
        %for monpolar
        S.Bipolar='No';
        S.channel=unique(S.bipole);

end



% Read the electrode names
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
        if ~strcmp(header.ChannelNames,'')
            for i1=1:length(S.channel)
                D.channels.label{i1,1}=header.ChannelNames{S.channel(i1)};
            end
        else
            tmp = spm_select(1, '\.txt$', 'Select file for names of channels');
            fid=fopen(tmp);
            n=0;
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                n=n+1;
                D.channels.label{n,1}=tline;
            end
            fclose(fid);
        end
    case{'Yes'}
        try
            S.bipole;
        catch
            P=spm_str_manip(Fdata,'h');
            tmp = spm_select(1, '\.txt$', 'Select txt file for bipolar derivations', {}, P);
            S.bipole=load(tmp);
        end
        for i1=1:size(S.bipole,1)
            if S.bipole(i1,1)==S.bipole(i1,2)
                D.channels.label{i1,1}=header.ChannelNames{S.bipole(i1,1)};
            else
                D.channels.label{i1,1}=[header.ChannelNames{S.bipole(i1,1)} '-' header.ChannelNames{S.bipole(i1,2)}];
            end
        end
end
%normalise names (lower case and p)
switch S.SEEG
    case{'Yes'}
        D.channels.label=lower(D.channels.label);
%         for i1=1:size(D.channels.label,1)
%             D.channels.label{i1,1}(findstr(D.channels.label{i1,1},''''))='p';
%         end
end            




% try
%     CreateTemplate=S.CreateTemplate;
% catch
%     CreateTemplate = spm_input('Create template ','+1','Yes|No');
% end
% switch CreateTemplate
%     case{'No'}
%         try
%             Fchannels = S.Fchannels;
%         catch
%             Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
%         end
%     case{'Yes'}
% %         P=spm_str_manip(Fdata,'h');
% %         filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
% %         fid=fopen(filename,'r');
% %         n=0;
% %         while 1
% %             tmp=fgetl(fid);
% %             if ~ischar(tmp)
% %                 break
% %             else
% %                 n=n+1;
% %                 Cnames{n,1}=tmp;
% %             end
% %         end
% %         fclose(fid);
%         Cnames=D.channels.label;
% %         for i1=1:size(Cnames,1)
% %             Cnames{i1,1}(find(Cnames{i1,1}==''''))='p';
% %         end
%         P=spm_str_manip(Fdata,'h');
%         try
%             filename = S.filenamePos;
%         catch
%             filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
%         end
%         Position=load(filename);
%         
%         %match between MRI and SEEG file
%         switch S.SEEG
%             case{'Yes'}
%                 Position2=zeros(size(Cnames,1),3);
%                 try
%                     filename = S.filenameName;
%                 catch
%                     filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
%                 end
%                 clear Name
%                 fid=fopen(filename);
%                 for i1=1:size(Position,1)
%                     Name{i1}=fgetl(fid);
%                 end
%                 fclose(fid);
%                 for i1=1:size(Cnames,1)
%                     [tmp1,tmp2]=strtok(Cnames{i1,1},'-');
%                     clear I1 I2
%                     for i2=1:length(Name)
%                         if strcmp(Name{i2},tmp1)
%                             I1=i2;
%                             break
%                         end
%                     end
%                     for i2=1:length(Name)
%                         if strcmp(Name{i2},tmp2(2:end))
%                             I2=i2;
%                             break
%                         end
%                     end
%                     Position2(i1,:)=mean(Position([I1 I2],:));
%                 end
%                 Position=Position2;
%         end
%                     
%         Cpos2=Position';
%         Cpos=Position(:,1:2)';
%         Cpos=Cpos-min(Cpos(:));
%         Cpos=Cpos./max(Cpos(:));
% 
% %         for i1=1:size(Cnames,1)
% %             tmp=find(Cnames{i1,1}=='p');
% %             tmp=tmp([length(tmp)/2 length(tmp)]);
% %             Cnames{i1,1}(tmp)='''';
% %         end
%         for i1=1:size(Cnames,1)
%             Cnames{i1,1}(find(Cnames{i1,1}=='-'))='';
%         end
%         
%         try
%             Species=S.Atlas;
%         catch
%             Species = spm_input('Species ','+1','Human|Rat|Mouse');
%         end
%         switch Species
%             case{'Rat','Mouse'}
%                 Cpos2=10*Cpos2;
%         end
% 
% 
%         Files=what(fullfile(spm('dir'),'EEGTemplates'));
%         ok=1;
%         while ok
%             try
%                 Filename=S.MontageName;
%                 ok=0;
%             catch
%                 Filename=spm_input('Name of montage', '+1','s');
%                 ok=0;
%                 for i1=1:length(Files.mat)
%                     if strcmp(lower(Files.mat{i1}),lower([Filename '.mat']))
%                         spm_input('Specified file already exists', '+1','d')
%                         ok=1;
%                         break
%                     end
%                 end
%             end
%         end
%         Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
%         Rxy=1;
%         Nchannels=length(Cnames);
%         if str2num(version('-release'))>=14
%             save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')
%         else
%             save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species')
%         end
% 
%         Fchannels=Filename;
%         Csetup = load(Fchannels);
% end


for i1=1:length(D.channels.label)
    D.channels.label{i1}(find(D.channels.label{i1}=='-'))='';
    
    %remove zeros
    a=D.channels.label{i1};
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    if length(v_num)>1
        if v_num(2)==v_num(1)+1 && strcmp(a(v_num(1)),'0')
            D.channels.label{i1}=a(setdiff(1:length(a),v_num(1)));
        end
    end

end



% % There doesn't seem to be information in the Deltamed-file which channel(s) is the reference
% % so ask for it now
% try
%     S.reference;
% catch
%     S.reference = spm_input('Input reference channel name', '+1', 's');
% end
S.reference=[];

try
    S.coarse;
catch
%     S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
    spm_input(sprintf('Sampling rate of acquisition = %d Hz',header.Radc),'+1','d');
    tmp=str2num(spm_input('Sampling rate of conversion [Hz]', '+1', 's',header.Radc));
    S.coarse=max([ceil(header.Radc/tmp) 1]);
    spm_input(sprintf('Effective sampling rate of conversion = %d Hz',header.Radc/S.coarse),'+1','d');
end

spm('Pointer','Watch'); drawnow;

D.Atlas=S.Atlas;

% if strcmp(CreateTemplate,'No')
%     Csetup = load(Fchannels);
% %     if length(S.channel)~=Csetup.Nchannels&~isempty(S.channel)
% %         error('Number of read channels does not match the template')
% %     end
% end


F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');


% file size
Dn = dir(Fdata);
Sf = Dn.bytes;

% Read data
fp = fopen(Fdata, 'r');
Data=fread(fp,[header.Nchannels,header.Nsamples],'int16');
fclose(fp);

if isfield(S,'channel')
    Data=Data(S.channel,:);
elseif isfield(S,'bipole')
    Datatmp=zeros(size(S.bipole,1),size(Data,2));
    for i1=1:size(S.bipole,1)
        if S.bipole(i1,1)==S.bipole(i1,2)
            Datatmp(i1,:)=Data(S.bipole(i1,1),:);
        else
            Datatmp(i1,:)=Data(S.bipole(i1,1),:)-Data(S.bipole(i1,2),:);
        end
    end
    Data=Datatmp;
end


% Read data and time strings
D.descrip.date = header.date;
D.descrip.time = header.time;

% Read number of channels
Nchannels = size(Data,1);

D.timeOnset=0;

% Read AD rate
D.Fsample=header.Radc;

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

%ECG channel
isecg=[];
for i1=1:length(D.channels)
    switch D.channels(i1).label
        case{'ecg','ecg1','ecg2','ekg1','ekg2'}
            D.channels(i1).type='ECG';
            isecg(end+1)=i1;
    end
end
%reorder to put ECG channels at the end
reorder=[setdiff(1:length(S.channel),isecg) isecg];
D.channels=D.channels(reorder);
Data=Data(reorder,:);


% % Map name of channels to channel order specified in channel template file
% for i = 1:length(D.channels)
% 	index = [];
% 	for j = 1:Csetup.Nchannels
% 		if ~isempty(find(strcmpi(D.channels(i).label, Csetup.Cnames{j})))
% 			index = [index j];
% 		end
% 	end
%     
% 	if isempty(index)
% %         disp('Imported channels are:')
% %         disp(D.channels.name')
% %         disp('Template channels are:')
% %         disp(Csetup.Cnames)
% % 		error(sprintf('No channel named %s found in channel template file.', D.channels.name{i}));
%         warning(sprintf('No channel named %s found in channel template file or in header.', D.channels(i).label));
% 	else
% 		% take only the first found channel descriptor
% 		D.channelOrder(i) = index(1);
%     end
% end

% D.data.fnamedat = [S.SaveFile '.dat'];
% D.data.datatype = 'float32-le';

D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
nsampl=Nsamples;
nchan=D.Nchannels;

%Events
%D.trials.events=header.evt;
D.trials.onset = 1/D.Fsample;

% if ~isempty(Ec)
%     for i1=1:max(Ec)
%         D.trials(i1).events=find(Ec==i1);
%         D.trials(i1).onset=Etsec(find(Ec==i1));
%     end
% else
%     D.trials.label = 'Undefined';
%     D.trials.events = [];
%     D.trials.onset = 1/D.Fsample;
% end

%Continuous
D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
% physically initialise file
D.data(end,end) = 0;
offset = 1;
nblocksamples = size(Data,2);
D.data(:, offset:(offset+nblocksamples-1)) = full(Data);


%Electrodes
for i1=1:length(D.channels)
    D.sensors.eeg.label{1,i1}=D.channels(i1).label;
%     D.sensors.eeg.label(i1)=D.channels(i1).label;
end
D.sensors.eeg.unit='mm';
D.sensors.eeg.pnt=NaN*zeros(length(D.channels),3);
D.sensors.eeg.type='eeg';
D.sensors.chantype='eeg';


D.Atlas=S.Atlas;

%--------- Create meeg object
[DirOut,NameOut,~]=fileparts(FileOut);
D.fname = [NameOut '.mat'];
D.path= DirOut;
D = meeg(D);
try
    D = events(D, 1, evt);
end
save(D)

spm('Pointer','Arrow');


function header=readheaderbin(File,Channels);

fid=fopen(File,'r');
for i1=1:16
    ligne = fgets(fid);
end
header.date=ligne(6:end);
ligne = fgets(fid);
header.time=ligne(6:end);
ligne = fgets(fid);
header.Radc=str2num(ligne(10:end));
ligne = fgets(fid);
ligne = fgets(fid);
ligne = fgets(fid);
header.Nsamples=str2num(ligne(19:end));
ligne = fgets(fid);
ligne = fgets(fid);
header.Nchannels=str2num(ligne(14:end));
ligne = fgets(fid);
ligne=ligne(10:end);
n=0;
while ligne~=-1
    n=n+1;
    [header.ChannelNames{n},ligne]=strtok(ligne,',');
    header.ChannelNames{n}=deblank(header.ChannelNames{n});
end
ligne = fgets(fid);
ligne = fgets(fid);
ligne=ligne(6:end);
[header.units,ligne]=strtok(ligne,','); 
for i1=1:7
    ligne = fgets(fid);
end
% n=0;
% while 1
%     ligne = fgets(fid);
%     if ~ischar(ligne) || length(ligne)<5, break, end
%     n=n+1;
%     [evt(n).time,ligne]=strtok(ligne,',');
%     evt(n).time=str2num(evt(n).time)/header.Radc;
%     [evt(n).type,ligne]=strtok(ligne,' ');
%     evt(n).type=deblank(evt(n).type(2:end));
%     if ~strcmp(evt(n).type,'Sélection')
%         [ligne,evt(n).duration]=strtok(ligne,'=');
%         evt(n).duration=str2num(evt(n).duration(2:end));
%     else
%         evt(n).value=n;
%     end
% end
% header.evt=evt;
fclose(fid);
