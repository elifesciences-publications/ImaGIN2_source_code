function Dnew = ImaGIN_BipolarMontage(S)
% Perform a longitudinal bipolar montage

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
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David

try
    Filename = S.Fname;
catch
    Filename = spm_select(Inf, '\.mat$', 'Select data file');
end

try
    FileOut=S.FileOut;
catch
    S.FileOut = Filename;
end

try
    SaveChannels = S.Save;
catch
    SaveChannels = 0;
end

try
    bipole = load(S.Bipole);
    
catch
    bipole = [];
end


for i0=1:size(Filename,1)
    T=deblank(Filename(i0,:));
    D=spm_eeg_load(T);
    
    Sensors=sensors(D,'EEG');
    Name=Sensors.label;
    Position=Sensors.elecpos;
    BadChannelsMono=badchannels(D);
    
    if isempty(bipole)
        Cpos2full=[];
        Cnamesfull={};
        bipole=[];
        for i1=1:length(Name)-1
            tmp1=Name{i1}(1:min([findstr(Name{i1},'0') findstr(Name{i1},'1') findstr(Name{i1},'2') findstr(Name{i1},'3') findstr(Name{i1},'4') findstr(Name{i1},'5') findstr(Name{i1},'6') findstr(Name{i1},'7') findstr(Name{i1},'8') findstr(Name{i1},'9')])-1);
            tmp2=Name{i1+1}(1:min([findstr(Name{i1+1},'0') findstr(Name{i1+1},'1') findstr(Name{i1+1},'2') findstr(Name{i1+1},'3') findstr(Name{i1+1},'4') findstr(Name{i1+1},'5') findstr(Name{i1+1},'6') findstr(Name{i1+1},'7') findstr(Name{i1+1},'8') findstr(Name{i1+1},'9')])-1);
            tmp3=str2num(Name{i1}(min([findstr(Name{i1},'0') findstr(Name{i1},'1') findstr(Name{i1},'2') findstr(Name{i1},'3') findstr(Name{i1},'4') findstr(Name{i1},'5') findstr(Name{i1},'6') findstr(Name{i1},'7') findstr(Name{i1},'8') findstr(Name{i1},'9')]):end));
            tmp4=str2num(Name{i1+1}(min([findstr(Name{i1+1},'0') findstr(Name{i1+1},'1') findstr(Name{i1+1},'2') findstr(Name{i1+1},'3') findstr(Name{i1+1},'4') findstr(Name{i1+1},'5') findstr(Name{i1+1},'6') findstr(Name{i1+1},'7') findstr(Name{i1+1},'8') findstr(Name{i1+1},'9')]):end));
            if strcmp(tmp1,tmp2)&tmp3+1==tmp4
                Cpos2full(:,end+1)=mean(Position([i1 i1+1],:))';
                Cnamesfull{1,end+1}=[Name{i1+1} Name{i1}];
                bipole=[bipole [i1+1;i1]];
            end          
        end
        if strcmp(Name{end},'ecg')
            Cpos2full(:,end+1)=mean(Position(end,:),1)';
            Cnamesfull{1,end+1}=Name{end};
            bipole=[bipole [length(Name);length(Name)]];
        end
    else
        bipole=bipole';
        Cpos2full=zeros(3,size(bipole,2));
        Cnamesfull={};
        for i1=1:size(bipole,2)
            Cpos2full(:,i1)=mean(Position(bipole(:,i1),:),1)';
            Cnamesfull{1,end+1}=[Name{bipole(1,i1)} '-' Name{bipole(2,i1)}];
        end
    end

    
    Sensors.chanpos=Cpos2full';
    Sensors.elecpos=Cpos2full';
    Sensors.label=Cnamesfull;
    Sensors.chanunit=Sensors.chanunit(bipole(1,:));
    Sensors.chantype=Sensors.chantype(bipole(1,:));
    
    %New bad channels
    BadChannelsBip=[];
    for i1=1:length(BadChannelsMono)
        [a1,a2]=find(bipole==BadChannelsMono(i1));
        BadChannelsBip=[BadChannelsBip a2'];
    end
    BadChannelsBip=unique(BadChannelsBip);    
    
    
    
    Data=D(bipole(1,:),:,:)-D(bipole(2,:),:,:);
    if strcmp(Name{end},'ecg')
        Data(end,:)=D(bipole(1,end),:,:);
    end
    
    %Save as a newfile
    try
        newname=FileOut;
        Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
    catch
        newname=['b' D.fname];
        Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
    end
    Dnew(:,:,:)=Data;
    Dnew=sensors(Dnew,'EEG',Sensors);
    Dnew=chanlabels(Dnew,1:size(Data,1),Sensors.label);
    Dnew=chantype(Dnew,1:size(Data,1),chantype(D,bipole(1,:)));
    if ~isempty(BadChannelsBip)
        Dnew = badchannels(Dnew,BadChannelsBip,1); % add badchannel index in meeg object
    end
    
    save(Dnew);
    
    if SaveChannels
        [dir,name,~]=fileparts(newname);
        fid = fopen(fullfile(dir, ['recordings_' name '.txt']),'w');
        recording_output = (Dnew.sensors('eeg').label);
        fprintf(fid,'%s\n',recording_output{:});
        fclose(fid);
    end
    
    % Save the list of channels in the output file
    ImaGIN_save_log(fullfile(Dnew), 'Bipolar montage:', chanlabels(Dnew));
end

