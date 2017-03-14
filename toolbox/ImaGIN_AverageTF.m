function D = ImaGIN_AverageTF(S)

%Average TF in a given frequency band (from 3D to 2D data)
%Olivier David

% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_converteeg2mat(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be either 'CNT' or 'BDF'
%    Mname     - char matrix of input file name(s)
%    Fchannels - String containing name of channel template file
%_______________________________________________________________________
%
% spm_eeg_converteeg2mat is a user interface to convert EEG-files from their
% native format to SPM's data format. This function assembles some
% necessary information before branching to the format-specific conversion
% routines.
% The user has to specify, by either using struct S or the GUI, a 'channel
% template file' that contains information about the (approximate) spatial
% positions of the channels.
%
% Output: The converted data are written to files. The header
% structs, but not the data, are returned in D as a cell vector of structs.
%_______________________________________________________________________
%
% Additional formats can be added by (i) extending the code below in a
% straightforward fashion, (ii) providing a new channel template file and
% (iii) adding the actual conversion routine to the SPM-code.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_TF_images.m 299 2005-11-15 15:25:17Z james $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TF',0);
try
    D = S.D;
catch
    D = spm_select(inf, '\.mat$', 'Select EEG mat file');
end
P = spm_str_manip(D, 'H');

try
    if size(D,1)==1
        D = spm_eeg_load(D);
    else
        DD=D;
        clear D
        for i1=1:size(DD,1)
            D{i1}= spm_eeg_load(deblank(DD(i1,:)));
        end
    end
catch
    error(sprintf('Trouble reading file %s', D));
end

if iscell(D)
    try
        Method = S.Method;
    catch
        Ctype = {
            'Mean',...
            'Median',...
            'ITC'};
        str   = 'Type of averaging ';
        Sel   = spm_input(str, 2, 'm', Ctype);
        Method = Ctype{Sel};
    end
    DD=D;clear D
    switch Method
        case 'Mean'
%             data=0;
%             for i1=1:length(DD)
%                 data=data+double(DD{i1}(:,:,:));
%             end
%             data=data./length(DD);
            data=zeros([size(DD{1}) length(DD)]);
            for i1=1:length(DD)
               D =DD{i1};
                if length(size(data))==4
                    data(:,:,:,i1)=double(D(:,:,:));
                else
                    data(:,:,:,:,i1)=double(D(:,:,:,:));
                end
            end
            data=jg_nanmean(data,length(size(data)));           
        case 'ITC'
            data=0;
            for i1=1:length(DD)
                data=data+exp(1i*double(DD{i1}(:,:,:)));
            end
            data=abs(data./length(DD));
        case 'Median'
            data=zeros([size(DD{1}) length(DD)]);
            for i1=1:length(DD)
               D =DD{i1};
                if length(size(data))==4
                    data(:,:,:,i1)=double(D(:,:,:));
                else
                    data(:,:,:,:,i1)=double(D(:,:,:,:));
                end
            end
            data=nanmedian(data,length(size(data)));           
    end
    try
        Name=S.NewName;
    catch
        Name=spm_input('Name of new file', '+1', 's');
    end
    D=DD{1};
    % D=rmfield(D,'data');
    if isempty(Name)%Assume namefiles are numbered, have the same events
        for i1=1:length(D.fnamedat)
            if ~strcmp(DD{1}.fnamedat(i1),DD{2}.fnamedat(i1))
                fn=D.fnamedat;
                i2=find(fn(i1+1:end)=='_');
                fnamedat=[fn(1:i1-1) 'Mean_' fn(i1+i2(1)+1:end)];
                break
            end
        end
    else
        fnamedat=[Name '.dat'];
    end
    D=clone(D,fnamedat,[D.nchannels D.Nfrequencies size(D,3) 1]);
    D(:,:,:)=data;
    save(D);
else
    if isfield(D, 'Nfrequencies');
        try
            fmt = S.fmt;
        catch
            spm_input('average over ...', 1, 'd')
            Ctype = {
                'electrodes',...
                'events',...
                'frequency'};
            str   = 'Average over which dimension';
            Sel   = spm_input(str, 2, 'm', Ctype);
            fmt = Ctype{Sel};
        end

        switch fmt
            case {'electrodes'}
                Ctype = {
                    'Dat',...
                    'Img'};
                str   = 'Generate ';
                Sel   = spm_input(str, 2, 'm', Ctype);
                fmt2 = Ctype{Sel};

                switch fmt2
                    case {'Dat'}
                        data=(sum(D(1:round(D.nchannels/2),:,:,:),1)+sum(D(round(D.nchannels/2)+1:end,:,:,:),1))./D.nchannels;
                        D=clone(D,['e' D.fnamedat], [1 size(D,2) size(D,3)]);
                        D.tf.channels=1;
                        D(:,:,:)=data;
                        save(D);

                    case {'Img'}
                        try
                            D.electrodes_of_interest = S.thresholds.elecs;
                        catch
                            str = 'electrodes[s]';
                            Ypos = -1;

                            while 1
                                if Ypos == -1
                                    [D.electrodes_of_interest, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                                else
                                    D.electrodes_of_interest = spm_input(str, Ypos, 'r', [], [1 Inf]);
                                end


                                t=1:D.nchannels;
                                tmp=[];
                                for en=D.electrodes_of_interest;
                                    if isempty(find(t==en))
                                        tmp=[tmp,en];
                                    end
                                end
                                if isempty(tmp) break, end
                            end
                        end
                        try
                            D.Nregion = S.region_no;
                        catch
                            str = 'region number';
                            Ypos = -1;

                            while 1
                                if Ypos == -1
                                    [D.Nregion, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                                else
                                    D.Nregion = spm_input(str, Ypos, 'r', [], [1 Inf]);
                                end
                                if ~isempty(D.Nregion) break, end
                                str = 'No data';
                            end
                        end
                        %number the different types
                        Events=events(D);
                        if(~isempty(Events))
                            try
                                Events(1).type;
                            catch
                                Events=Events{1};
                            end
                        end
                        Types={};
                        for i1=1:length(Events)
                            trouve=0;
                            for i2=1:length(Types)
                                if ~strcmp(Types{i2},Events(i1).type) && i2==length(Types) && ~trouve
                                    Types{end+1}=Events(i1).type;
                                elseif strcmp(Types{i2},Events(i1).type)
                                    trouve=1;
                                end
                            end
                            if isempty(Types)
                                Types{1}=Events(1).type;
                            end
                        end
                        
                        for i = 1 : length(Types)
                            Itrials = find(Events.type == Events(i).type);
                            cd(D.path)
                            dname = sprintf('%dROI_TF_trialtype%d', D.Nregion, Events(i).types);
                            [m, sta] = mkdir(dname);
                            cd(dname);

                            for l = Itrials
                                % if single trial data make new directory for single trials,
                                % otherwise just write images to trialtype directory
                                if D.ntrials~= Types
                                    % single trial data
                                    dname = sprintf('trial%d.img', l);
                                    fname = dname;
                                    [m, sta] = mkdir(dname);
                                    cd(dname);
                                else
                                    fname = 'average.img';
                                end
                                data=squeeze(mean(D(D.electrodes_of_interest,:,:,i),1));
                                V.fname = fname;
                                V.dim = [D.Nfrequencies D.nsamples  1 ];
                                V.dt=[spm_type('float64') 0]; %%%check later with john
                                V.mat = eye(4);
                                V.pinfo = [1 0 0]';
                                spm_write_vol(V, data); % d is data
                            end
                        end
                end

            case {'frequency'}
                try
                    D.Frequency_window = S.freqs;
                    Ypos = -1;
                    while 1
                        if Ypos == -1
                            Ypos = '+1';
                        end
                        inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2));
                        if ~isempty(inds) break, end
                        str = 'No data in range';
                    end
                catch
                    str = 'Frequency window';

                    Ypos = -1;
                    while 1
                        if Ypos == -1
                            Ypos = '+1';
                        end
                        [D.Frequency_window, Ypos] = spm_input(str, Ypos, 'r', [], 2);

                        inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2))
                        if ~isempty(inds) break, end
                        str = 'No data in range';
                    end
                end
                data=squeeze(mean(D(:,inds,:,:),2));
                D=clone(D,['F' num2str(D.Frequency_window(1)) '_' num2str(D.Frequency_window(2)) '_' D.fnamedat], [size(D,1) size(D,3) 1]);
                D(:,:)=data;
                %          D(:,3)=[];
                D=rmfield(D,'Nfrequencies');
                %      probl          if isfield(D.tf,'time')
                %                     D=time(D,D.tf.time);
                %                 end
                D=rmfield(D,'tf');
                save(D);

            case {'events'}
                %number the different types
                Events=events(D);
                if(~isempty(Events))
                    try
                        Events(1).type;
                    catch
                        Events=Events{1};
                    end
                end
                Types={};
                for i1=1:length(Events)
                    trouve=0;
                    for i2=1:length(Types)
                        if ~strcmp(Types{i2},Events(i1).type) && i2==length(Types) && ~trouve
                            Types{end+1}=Events(i1).type;
                        elseif strcmp(Types{i2},Events(i1).type)
                            trouve=1;
                        end
                    end
                    if isempty(Types)
                        Types{1}=Events(1).type;
                    end
                end
                data=squeeze(mean(D(:,:,:,:),4));
%                Events.time=Events(1).time;
%                 value=Events(i1).value;
%                 D.events=rmfield(Events,'value');
%                 Events(1).value=value;
%                 Event.type='1';
                D=clone(D,['m_' D.fnamedat], [D.nchannels D.Nfrequencies D.nsamples 1]);
                D=events(D,1,Events);
                D(:,:,:)=data;
                save(D);
        end
    else
        error('No time frequency data');
    end
end
