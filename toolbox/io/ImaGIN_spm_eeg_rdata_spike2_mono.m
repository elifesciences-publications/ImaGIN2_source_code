function D = ImaGIN_spm_eeg_rdata_spike2_mono(S)
% converts EEG data from Spike2- to SPM-format
% FORMAT D = spm_eeg_rdata(S)
%
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of Spike2-file
% Fchannels   - filename of channel template
% reference   - name of reference channel(s)
% channel     - index of channel(s) to read
%_______________________________________________________________________
%
% ImaGIN_spm_eeg_rdata_spike2_mono reads a continuous *.smr Spike2 file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rdata.m 317 2005-11-28 18:31:24Z stefan $

Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Finter)


try
    Fdata = S.Fdata;
catch
    Fdata = spm_select(1, '\.smr$', 'Select Spike2 file');
end

fp = fopen(Fdata, 'r');
%Read header
Head=SONFileHeader(fp);


try
    FileOut=S.FileOut;
catch
    FileOut = spm_input('Name of converted file', '+1', 's');
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

% try
%     CreateTemplate=S.CreateTemplate;
% catch
%     CreateTemplate = spm_input('Create template ','+1','Yes|No');
% end
% 
% switch CreateTemplate
%     case{'No'}
%         try
%             Fchannels = S.Fchannels;
%         catch
%             Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
%         end
% end
% 


% % There doesn't seem to be information in the CNT-file which channel(s) is the reference
% % so ask for it now
% try
%     S.reference;
% catch
%     S.reference = spm_input('Input reference channel name', '+1', 's');
% end
S.reference=[];



% try
%     Bipolar=S.Bipolar;
% catch
%     Bipolar = spm_input('Bipolar derivations? ','+1','Yes|No');
% end
% switch Bipolar
%     case{'No'}
        try
            S.channel;
        catch
            S.channel = str2num(spm_input('Index of channels to read', '+1', 's'));
        end
%     case{'Yes'}
%         try
%             S.bipole=load(S.bipole);
%         catch
%             P=spm_str_manip(Fdata,'h');
%             tmp = spm_select(1, '\.txt$', 'Select txt file for bipolar derivations', {}, P);
%             S.bipole=load(tmp);
%         end
%         S.channel=unique(S.bipole);
%         %         for i1=1:size(S.bipole,1)
%         %             D.channels.name{i1}=[header.ChannelNames{S.bipole(i1,1)} '-' header.ChannelNames{S.bipole(i1,2)}];
%         %         end
% end

try
    S.epochlength;
catch
    S.epochlength = str2num(spm_input('Max epoch length [s] ', '+1', 's'));
end

% try
%     Continuous=S.continuous;
% catch
%     Continuous = spm_input('Continuous recordings? ','+1','Yes|No');
%     if strcmp(Continuous,'Yes')
%         Continuous=1;
%     else
%         Continuous=0;
%     end
% end

try
    S.coarse;
catch
    S.coarse = str2num(spm_input('Reading downsampling factor ', '+1', 's',1));
end

% try
%     S.pre;
% catch
%     S.pre = spm_input('Prefixe filename ', '+1', 's');
% end




spm('Pointer','Watch'); drawnow;

% Read data
if isempty(S.epochlength)
    [Data,Param.ChannelNumber,Evt]=ReadSMR(Fdata,S.channel,S.coarse);

    D = [];

    % Read the electrode names
%     switch Bipolar
%         case{'No'}
            D.channels.label = Data.name;
%         case{'Yes'}
%             Datatmp=zeros(size(S.bipole,1),size(Data.data,2));
%             for i1=1:size(S.bipole,1)
%                 if S.bipole(i1,1)==S.bipole(i1,2)
%                     D.channels.label{i1}=Data.name{find(S.channel==S.bipole(i1,1))};
%                     Datatmp(i1,:)=Data.data(find(S.channel==S.bipole(i1,1)),:);
%                 else
%                     D.channels.label{i1}=[Data.name{find(S.channel==S.bipole(i1,1))} '-' Data.name{find(S.channel==S.bipole(i1,2))}];
%                     Datatmp(i1,:)=Data.data(find(S.channel==S.bipole(i1,1)),:)-Data.data(find(S.channel==S.bipole(i1,2)),:);
%                 end
%             end
%             Data.data=Datatmp;
%             clear Datatmp
%     end

    %spikes
    if isfield(Data,'spike')
        D.spike=Data.spike;
        %     D.channels.name{end+1:end+length(Data.spike.name)}=Data.spike.name{:};
    end

    
%             if strcmp(CreateTemplate,'Yes')
%             Cnames=D.channels.name';
%             filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
%             Position=load(filename);
%             if size(Position,1)~=length(Cnames)
%                 error('Number of electrode positions and of electrode names is not the same')
%             end
% 
%             Cpos2=Position';
%             Cpos=Position(:,1:2)';
%             Cpos=Cpos-min(Cpos(:));
%             Cpos=Cpos./max(Cpos(:));
% 
%             Species = S.Atlas;
%             switch Species
%                 case{'Rat','Mouse'}
%                     Cpos2=10*Cpos2;
%             end
% 
%             Files=what(fullfile(spm('dir'),'EEGTemplates'));
%             ok=1;
%             while ok
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
%             Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
%             Rxy=1;
%             Nchannels=length(Cnames);
%             save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')
% 
%             Fchannels=Filename;
%         end
        
    F = spm_str_manip(Fdata, 'rt');
    P = spm_str_manip(Fdata, 'H');
    

    % file size
    Dn = dir(Fdata);
    Sf = Dn.bytes;





%     Csetup = load(Fchannels);


    % Read data and time strings
    D.descrip.date = [num2str(Head.timeDate.Year) '/' num2str(Head.timeDate.Detail(6)) '/' num2str(Head.timeDate.Detail(5))];
    D.descrip.time = [num2str(Head.timeDate.Detail(4)) '/' num2str(Head.timeDate.Detail(3)) '/' num2str(Head.timeDate.Detail(2))];

    % Read number of channels
    Nchannels = size(Data.data,1);

    % Read AD rate
    D.Fsample = inv(Data.time(2)-Data.time(1));

%     %Events
%     Nevent = str2num(spm_input('Number of event types', '+1', 's',0));
%     Nevents=0;
    Ec=[];
%     Et=[];
%     for i1=1:Nevent
%         try
%             D.trials.label{i1}=S.event(i1); %I changed S.event{i1}
%         catch
%             D.trials.label{i1}=spm_input(sprintf('Name of event %d',i1), '+1', 's');
%         end
% 
%         try
%             tmp=deblank(S.event_file(i1,:));
%         catch
%             tmp = spm_select(1, '\.txt$', sprintf('Select txt file with the timing (sec) of event %d',i1));
%         end
%         tmp=load(tmp);
%         Nevents=Nevents+length(tmp);
%         Ec=[Ec;i1*ones(length(tmp),1)];
%         for i2=1:length(tmp)
%             Et=[Et;unique(find(abs(Data.time-tmp(i2))==min(abs(Data.time-tmp(i2)))))];
%         end
%     end
%     Etsec=Data.time(Et);

    % Number of expected samples
    Nsamples = length(Data.time);

    % store only relative name of channel template file
    D.channels = repmat(struct('bad', 0,'label',D.channels.label), 1,1);

    %Time zero (OD)
    [tmp D.TimeZero]=min(abs(Data.time));
    

    D.time=Data.time;

    D.Atlas=S.Atlas;

%     % Map name of channels to channel order specified in channel template file
%     for i = 1:length(D.channels)
%         index = [];
%         for j = 1:Csetup.Nchannels
%             if ~isempty(find(strcmpi(D.channels(i).label, Csetup.Cnames{j})))
%                 index = [index j];
%             end
%         end
% 
%         if isempty(index)
%             warning(sprintf('No channel named %s found in channel template file.', D.channels(i).label));
%         else
%             % take only the first found channel descriptor
%             D.channels(i).order = index(1);
%         end
%     end


        
        D.Nchannels = Nchannels;
        D.Nsamples = Nsamples;
        nsampl=Nsamples;
        nchan=Nchannels;


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
        nblocksamples = size(Data.data,2);
        D.data(:, offset:(offset+nblocksamples-1)) = full(Data.data);


%         %Electrodes
%         try
%             D.sensors.eeg.pnt=Csetup.Cpos2(:,ChannelOrder)';
%         catch
%             D.sensors.eeg.pnt=Csetup.Cpos2(:,:)';
%         end
%         for i1=1:length(D.channels)
%             D.sensors.eeg.label{1,i1}=D.channels(i1).label;
%         end
%         D.sensors.eeg.unit='mm';
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
        D.fname = [FileOut '.mat'];

        D = meeg(D);

        if ~isempty(Evt)
            D = events(D, 1, Evt);
        end
        %     D=timeonset(D,D.time(D.TimeZero));

        %spikes
        if isfield(Data,'spike')
            D.spike=Data.spike;
            %     D.channels.name{end+1:end+length(Data.spike.name)}=Data.spike.name{:};
        end
        save(D)

        spm('Pointer','Arrow');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    fid=fopen(Fdata);
    Info=SONChannelInfo(fid,S.channel(1));
    fclose(fid);
    time=(0:Info.maxChanTime-1)./Info.idealRate;
    NBlocks=ceil(max(time)/S.epochlength);
    ReadTime.NSamples=ceil(S.epochlength*Info.idealRate);
    for i0=1:NBlocks
        ReadTime.Start=(i0-1)*ReadTime.NSamples+1;
        ReadTime.Stop=i0*ReadTime.NSamples;
        [Data,Param.ChannelNumber]=ImaGIN_ReadSMR(Fdata,S.channel,S.coarse,ReadTime);

        D = [];

        % Read the electrode names
        switch Bipolar
            case{'No'}
                D.channels.label = Data.name;
            case{'Yes'}
                Datatmp=zeros(size(S.bipole,1),size(Data.data,2));
                for i1=1:size(S.bipole,1)
                    if S.bipole(i1,1)==S.bipole(i1,2)
                        D.channels.label{i1}=Data.name{find(S.channel==S.bipole(i1,1))};
                        Datatmp(i1,:)=Data.data(find(S.channel==S.bipole(i1,1)),:);
                    else
                        D.channels.label{i1}=[Data.name{find(S.channel==S.bipole(i1,1))} '-' Data.name{find(S.channel==S.bipole(i1,2))}];
                        Datatmp(i1,:)=Data.data(find(S.channel==S.bipole(i1,1)),:)-Data.data(find(S.channel==S.bipole(i1,2)),:);
                    end
                end
                Data.data=Datatmp;
                clear Datatmp
        end

        

        
        if strcmp(CreateTemplate,'Yes') && i0==1
            Cnames=D.channels.name';
            filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
            Position=load(filename);
            if size(Position,1)~=length(Cnames)
                error('Number of electrode positions and of electrode names is not the same')
            end

            Cpos2=Position';
            Cpos=Position(:,1:2)';
            Cpos=Cpos-min(Cpos(:));
            Cpos=Cpos./max(Cpos(:));

            Species = S.Atlas;
            switch Species
                case{'Rat','Mouse'}
                    Cpos2=10*Cpos2;
            end

            Files=what(fullfile(spm('dir'),'EEGTemplates'));
            ok=1;
            while ok
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
            Filename=fullfile(spm('dir'), 'EEGtemplates',Filename);
            Rxy=1;
            Nchannels=length(Cnames);
            save(Filename,'Cnames','Cpos','Cpos2','Nchannels','Rxy','Species','-V6')

            Fchannels=Filename;
        end

        
        F = spm_str_manip(Fdata, 'rt');
        P = spm_str_manip(Fdata, 'H');


        % file size
        Dn = dir(Fdata);
        Sf = Dn.bytes;



        Csetup = load(Fchannels);
        % if length(D.channels.name)~=Csetup.Nchannels
        %     error('Number of read channels does not match the template')
        % end


        % Read data and time strings
        D.descrip.date = [num2str(Head.timeDate.Year) '/' num2str(Head.timeDate.Detail(6)) '/' num2str(Head.timeDate.Detail(5))];
        D.descrip.time = [num2str(Head.timeDate.Detail(4)) '/' num2str(Head.timeDate.Detail(3)) '/' num2str(Head.timeDate.Detail(2))];

        % Read number of channels
        Nchannels = size(Data.data,1);

        % Read AD rate
        D.Fsample = inv(Data.time(2)-Data.time(1));

        %Events
%         Nevent = str2num(spm_input('Number of event types', '+1', 's',0));
        Nevent = 0;
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
                Et=[Et;unique(find(abs(Data.time-tmp(i2))==min(abs(Data.time-tmp(i2)))))];
            end
        end
        Etsec=Data.time(Et);

        % Number of expected samples
        Nsamples = length(Data.time);

        % store only relative name of channel template file
        D.channels = repmat(struct('bad', 0,'label',D.channels.label), 1, 1);

        %Time zero (OD)
        [tmp D.TimeZero]=min(abs(Data.time));

        D.time=Data.time;


        D.Atlas=S.Atlas;

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
        
        D.Nchannels = Nchannels;
        D.Nsamples = Nsamples;
        nsampl=Nsamples;
        nchan=Nchannels;


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
        D.data = file_array([S.NameOut '.dat'], [nchan nsampl], 'float32-le');
        % physically initialise file
        D.data(end,end) = 0;
        offset = 1;
        nblocksamples = size(Data.data,2);
        D.data(:, offset:(offset+nblocksamples-1)) = full(Data.data);


%         %Electrodes
%         try
%             D.sensors.eeg.pnt=Csetup.Cpos2(:,ChannelOrder)';
%         catch
%             D.sensors.eeg.pnt=Csetup.Cpos2(:,:)';
%         end
%         for i1=1:length(D.channels)
%             D.sensors.eeg.label{1,i1}=D.channels(i1).label;
%         end
%         D.sensors.eeg.unit='mm';
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

        D = meeg(D);
        D = path(DirOut);
        
        %spikes
        if isfield(Data,'spike')
            D.spike=Data.spike;
            %     D.channels.name{end+1:end+length(Data.spike.name)}=Data.spike.name{:};
        end
        save(D)
    end
end
 spm('Pointer','Arrow');


