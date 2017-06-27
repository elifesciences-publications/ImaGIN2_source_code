function ImaGIN_Crop(P)
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

try
    Job = P.Job;
catch
    Job = spm_input('Method ','1','Epoch|Manual');
end

switch Job
    case{'Epoch'}
        try
            NewFile = P.NewFile;
        catch
            NewFile = spm_input('Create a file per event ','1','Yes|No',[1 0]);
        end
        try
            t=P.Fname;
        catch
            t = spm_select(1, '\.mat$', 'Select data file');
        end
        D=spm_eeg_load(t);
        D.NewFile=NewFile;
        try
            D.EventRef=P.EventRef;
        end
        try
            D.OffsetStart=P.OffsetStart;
        end
        try
            D.OffsetEnd=P.OffsetEnd;
        end
        ImaGIN_CropEpoch(D);

    case{'Manual'}
        try
            t = P.Fname;
        catch
            t = spm_select(1, '\.mat$', 'Select data file');
        end
        D = spm_eeg_load(t);
        try
            ImaGIN_CropManual(D,t,P);
        catch
            ImaGIN_CropManual(D,t);
        end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_CropManual(D,t,S)
    try
        EventRefName=S.EventStart;
    catch
        allEvt = unique({D.events.type});
        disp(['Available event types: ', sprintf('%s ', allEvt{:})]);
        EventRefName=spm_input('Type of start event', '+1', 's');
    end
    if isempty(EventRefName)
        EventStart=0;
    end
    Events=D.events;
    try
        EventStart=Events(S.numbStart).time-timeonset(D);
    catch
        n=1;
        for i=1:length(Events)
            if strcmp(EventRefName,strvcat(Events(i).type))
                EventStart(n)=Events(i).time-timeonset(D);
                n=n+1;
            end
        end
    end
    EventStart=min(ceil(EventStart*D.fsample));
    try
        OffsetStart=S.OffsetStart;
    catch
        OffsetStart=spm_input('Duration [sec] before start event', '+1', 'r',0);
    end
    OffsetStart=round(OffsetStart*D.fsample);
    Time(1)=max([1 EventStart-OffsetStart]);

    try
        EventRefNameEnd=S.EventEnd;
    catch
        EventRefNameEnd=spm_input('Type of end event', '+1', 's');
    end
    if isempty(EventRefNameEnd)
        EventEnd=0;
    end
    try
        EventEnd=Events(S.numbEnd).time-timeonset(D);
    catch
        n=1;
        for i=1:length(Events)
            if strcmp(EventRefNameEnd,strvcat(Events(i).type))
                EventEnd(n)=Events(i).time-timeonset(D);
                n=n+1;
            end
        end
    end
    EventEnd=max(ceil(EventEnd*D.fsample));
    try
        OffsetEnd=S.OffsetEnd;
    catch
        OffsetEnd=spm_input('Duration [sec] after end event', '+1', 'r',0);
    end
    OffsetEnd=round(OffsetEnd*D.fsample);
    Time(2)=min([D.nsamples EventEnd+OffsetEnd]);

    %Crop
    Index=Time(1):Time(2);
    if isfield(D,'time')
        D.time=D.time(Index);
    end
    d = D(:, Index, 1);

    %Save as a newfile
    try
        Dnew = clone(D, S.FileOut, [D.nchannels length(Index), 1]);
    catch
        try
            Prefix=S.Prefix;
        catch
            Prefix=spm_input('Prefix of new file', '+1', 's');
        end
        Dnew = clone(D, [Prefix fname(D)], [D.nchannels length(Index), 1]);
    end

    Dnew(:, :, 1) = d;
    Dnew = events(Dnew, 1, select_events(Events,[Time(1)/D.fsample+timeonset(D)  Time(2)/D.fsample+timeonset(D)]));
    Dnew = timeonset(Dnew, Time(1)./D.fsample+D.timeonset);

    save(Dnew);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_CropEpoch(D)
    try
        EventRefName=D.EventRef;
    catch
        allEvt = unique({D.events.type});
        disp(['Available event types: ', sprintf('%s ', allEvt{:})]);
        EventRefName=spm_input('Type of reference event', '+1','s');
    end
    Events=D.events;
    n=1;
    for i=1:length(Events)
        if strcmp(EventRefName,strvcat(Events(i).type))
            EventRef(n)=Events(i).time;
            n=n+1;
        end
    end
    EventRef=indsample(D,EventRef);

    try
        OffsetStart=D.OffsetStart;
    catch
        OffsetStart=spm_input('Duration [sec] before reference event', '+1', 'r');
    end
    OffsetStartName=OffsetStart;
    OffsetStart=ceil(OffsetStart*D.fsample);
    try
        OffsetEnd=D.OffsetEnd;
    catch
        OffsetEnd=spm_input('Duration [sec] after reference event', '+1', 'r');
    end
    OffsetEndName=OffsetEnd;
    OffsetEnd=ceil(OffsetEnd*D.fsample);

    n=0;
    DD=D;
    try
        time=DD.time;
    catch
        time = 0 : (1/DD.fsample) : ((DD.nsamples-1)/DD.fsample);
    end

    if D.NewFile
        for i1=1:length(EventRef)
            if (EventRef(i1)-OffsetStart>=1) && (EventRef(i1)-OffsetStart<=DD.nsamples) && (EventRef(i1)+OffsetEnd>=1) && (EventRef(i1)+OffsetEnd<=DD.nsamples)
                n=n+1;
                D=DD;
                Time(1)=EventRef(i1)-OffsetStart;
                Time(2)=EventRef(i1)+OffsetEnd;
                Index=Time(1):Time(2);

                %Save as a newfile
                if n>99
                    Prefix=[EventRefName '_' sprintf('%d-%d_%d',round(OffsetStartName),round(OffsetEndName),n)];
                elseif n>9
                    Prefix=[EventRefName '_' sprintf('%d-%d_0%d',round(OffsetStartName),round(OffsetEndName),n)];
                else
                    Prefix=[EventRefName '_' sprintf('%d-%d_00%d',round(OffsetStartName),round(OffsetEndName),n)];
                end

                Dnew = clone(D, [Prefix '_' fname(D)], [D.nchannels length(Index), 1]);
                d = D(:, Index, 1); 
                Dnew(:, :, 1) = d;
                Dnew = events(Dnew, 1, select_events(D.events,[Time(1)/D.fsample+time(1)  Time(2)/D.fsample+time(1)]));
                Dnew = timeonset(Dnew, D.time(EventRef(i1)-OffsetStart)-D.time(EventRef(i1)));

                Events=events(Dnew);
                for i0=1:length(Events)
                    Events(i0).time = Events(i0).time-D.time(EventRef(i1));
                end
                if isfield(D,'spike')
                    for i2=1:length(D.spike.timings)
                        D.spike.timings{i2}=D.spike.timings{i2}-D.time(EventRef(i2));
                    end
                end
                Dnew = events(Dnew, 1, Events);

                save(Dnew);
            end
        end

    else
        D=DD;
        %Save as a newfile
        Prefix=[EventRefName '_' sprintf('%d-%d',round(OffsetStartName),round(OffsetEndName))];
        ntrials=length(EventRef);
        Dnew = clone(D, [Prefix '_' fname(D)], [D.nchannels OffsetStart+OffsetEnd+1, ntrials]);
        for i1=1:ntrials
            if (EventRef(i1)-OffsetStart>=1) && (EventRef(i1)-OffsetStart<=DD.nsamples) && (EventRef(i1)+OffsetEnd>=1) && (EventRef(i1)+OffsetEnd<=DD.nsamples)
                n=n+1;
                Time(1)=EventRef(i1)-OffsetStart;
                Time(2)=EventRef(i1)+OffsetEnd;
                Index=Time(1):Time(2);
                d = D(:, Index, 1); 
                Dnew(:, :, i1) = d;
                Dnew = events(Dnew, i1, select_events(D.events,[Time(1)/D.fsample  Time(2)/D.fsample]));
                Dnew = timeonset(Dnew, Time(1)./D.fsample+D.timeonset);
                if isfield(DD,'spike')
                    for i2=1:length(DD.spike.timings)
                        tmp=find(DD.spike.timings{i2}>=time(EventRef(i1))-OffsetStart/D.fsample&DD.spike.timings{i2}<=time(EventRef(i1))+OffsetEnd/D.fsample);
                        D.spike.timings{i2,n}=DD.spike.timings{i2}(tmp)-time(EventRef(i1));
                        D.spike.markers{i2,n}=DD.spike.markers{i2}(tmp,:);
                        D.spike.adc{i2,n}=DD.spike.adc{i2}(tmp,:);
                    end
                end
                save(Dnew);
            else
                error('Duration before/after event is too long');
            end
        end
    end
end


%==========================================================================
% Utility function to select events according to time segment
function event = select_events(event, timeseg)
    if ~isempty(event)
        [time ind] = sort([event(:).time]);
        selectind = ind(time >= timeseg(1) & time <= timeseg(2));
        event = event(selectind);
    end
end


