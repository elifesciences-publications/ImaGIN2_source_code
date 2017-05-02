function D = ImaGIN_spm_eeg_rdata_elan(S)
% converts EEG data from Spike2- to SPM-format
% FORMAT D = spm_eeg_rdata(S)
%
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of Spike2-file
% Fchannels   - filename of channel template
% reference   - name of reference channel(s)
% channel     - index of channel(s) to read

% Stefan Kiebel
% $Id: spm_eeg_rdata.m 317 2005-11-28 18:31:24Z stefan $

Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Finter)


try
    Fdata = S.Fdata;
catch
    Fdata = spm_select(1, '\.eeg$', 'Select Spike2 file');
end

try
    FileOut=S.FileOut;
catch
    FileOut=spm_input('Name of converted file', '+1', 's');
end

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
            S.bipole=load(S.bipole);
        catch
            P=spm_str_manip(Fdata,'h');
            tmp = spm_select(1, '\.txt$', 'Select txt file for bipolar derivations', {}, P);
            S.bipole=load(tmp);
        end
        S.channel=unique(S.bipole);
        %         for i1=1:size(S.bipole,1)
        %             D.channels.name{i1}=[header.ChannelNames{S.bipole(i1,1)} '-' header.ChannelNames{S.bipole(i1,2)}];
        %         end
end


S.Atlas='Human';

spm('Pointer','Watch'); drawnow;

% Read data
    
[m_data,m_events,v_label_selected,s_fs,s_nb_samples_all,s_nb_channel_all,v_label_all,v_channel_type_all,v_channel_unit_all,str_ori_file1,str_ori_file2] = eeg2mat(Fdata,1,'all','all');
Data.data=m_data;
Data.name=v_label_selected;
S.channel=1:8;
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
CreateTemplate='No';
    
    if strcmp(CreateTemplate,'Yes')
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
    
    try
        S.NameOut;
    catch
        S.NameOut = spm_input('Name of converted file', '+1', 's');
    end
  

    % file size
    Dn = dir(Fdata);
    Sf = Dn.bytes;





%     Csetup = load(Fchannels);

    % Read number of channels
    Nchannels = size(Data.data,1);

    % Read AD rate
    D.Fsample = s_fs;

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
    Nsamples = s_nb_samples_all;

    % store only relative name of channel template file
    D.channels = repmat(struct('bad', 0,'label',D.channels.label), 1,1);    

    D.time=[0:Nsamples-1]./D.Fsample;

    D.Atlas=S.Atlas;


%        D.data.fnamedat = [Fdata(1:end-4) '.dat'];
%         D.data.datatype = 'float32-le';
        
        D.Nchannels = Nchannels;
        D.Nsamples = Nsamples;
        nsampl=Nsamples;
        nchan=Nchannels;


            D.trials.label = 'Undefined';
            D.trials.events = [];
            D.trials.onset = 1/D.Fsample;

        %Continuous
        D.data = file_array([FileOut '.dat'], [nchan nsampl], 'float32-le');
        % physically initialise file
        D.data(end,end) = 0;
        offset = 1;
        nblocksamples = size(Data.data,2);
        D.data(:, offset:(offset+nblocksamples-1)) = full(Data.data);


        %Electrodes
        try
            D.sensors.eeg.pnt=zeros(nchan,3);
        catch
            D.sensors.eeg.pnt=Csetup.Cpos2(:,:)';
        end
        for i1=1:length(D.channels)
            D.sensors.eeg.label{1,i1}=D.channels(i1).label;
        end
        D.sensors.eeg.unit='mm';


        D.Atlas=S.Atlas;

        %--------- Create meeg object
        D.fname = [Fdata(1:end-4) '.mat'];
        [DirOut,~,~]=fileparts(FileOut);
        D.path = DirOut;
        D = meeg(D);

        if ~isempty(m_events)
        evt=[];
%                 D.trials(n).onset=m_events(tmp,1)/D.Fsample; %in sec
                for i1=1:size(m_events,1)
                    evt(i1).type  = num2str(m_events(i1,2));
                evt(i1).time = m_events(i1,1)/D.fsample;
                evt(i1).value= m_events(i1,2);
                end
            D = events(D, 1, evt);
        end
            

        save(D)


        spm('Pointer','Arrow');


