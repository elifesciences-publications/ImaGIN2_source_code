function ImaGIN_Events(S)
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
    t = S.Fname;
catch
    t = spm_select(1, '\.mat$', 'Select data file');
end
if isempty(t)
    return;
end

try
    Action=S.Action;
catch
    Action = spm_input('Events ',1,'Add|Remove');
end

if strcmp(Action,'Add')
    try
        ImaGIN_EventsAdd(t,S);
    catch
        ImaGIN_EventsAdd(t);
    end
elseif strcmp(Action,'Remove')
    try
        ImaGIN_EventsRemove(t,S);
    catch
        ImaGIN_EventsRemove(t);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_EventsAdd(Filename,S)
    try
        NeventNew=S.Nevent;
    catch
        NeventNew=spm_input('Number of events to add', '+1', 'r',1);
    end

    for i1=1:NeventNew
        try
            NewName{i1}=S.EventName{i1};
        catch
            NewName{i1}=spm_input(sprintf('Type of event %d',i1), '+1', 's');
        end
        if isfield(S, 'Timing') && ~isempty(S.Timing)
            Timing{i1} = S.Timing{i1};
        else
            try
                tmp = S.EventFileName{i1};
            catch
                tmp = spm_select(1, '\.txt$', sprintf('Select txt file with the timing (sec) of event %d',i1));
            end
            Timing{i1}=load(tmp);
        end
    end

    try
        FileOut = S.FileOut;
    catch
        FileOut = Filename;
    end

    D=ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew);

    Direc=spm_str_manip(Filename,'h');
    File=spm_str_manip(Filename,'t');
    E=what(Direc);
    ok=0;

    D2=clone(D,FileOut, [D.nchannels D.nsamples D.ntrials]);
    D2(:,:,:)=D(:,:,:);
    save(D2);

    if ~isempty(E)
        for i1=1:length(E.mat)
            if strcmp(E.mat{i1},['t1_' File])
                ok=1;
                break
            end
        end
    end

    if ok
        Filename=fullfile(Direc,['t1_' File]);
        ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew);
        Filename=fullfile(Direc,['t2_' File]);
        ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew);
        Filename=fullfile(Direc,['t1int_' File]);
        ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew);
        Filename=fullfile(Direc,['t2int_' File]);
        ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=ImaGIN_EventsAdd_Subfunction(Filename,NewName,Timing,NeventNew)
    D=spm_eeg_load(Filename);

    NeventOld = size(D.events,2);
    if isempty(events(D))
        NeventOld=0;
    end

    evt = D.events;
    if ~isfield(evt,'type')
        NeventOld=0;
        clear evt
    end

    n=NeventOld;
    for i0=1:NeventNew
        for i1=1:length(Timing{i0})
            n=n+1;
            evt(n).type  = NewName{i0};
            evt(n).time  = Timing{i0}(i1);
            evt(n).value = NeventOld+i0;
        end
        % Add log entry
        ImaGIN_save_log(fullfile(D), sprintf('Added %dx event %s', length(Timing{i0}), NewName{i0}));
    end
    % This assigns these events to the first trials (the only one if you have continuous data)
    D = events(D, 1, evt);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_EventsRemove(Filename,S)
    D = spm_eeg_load(Filename);

    try
        EventRemove = S.EventName;
    catch
        EventRemove = spm_input('Type of events to remove (*=ALL)', '+1', 's');
    end

    % Get events in file
    Event = D.events;
    nEvent = size(Event,2);
    % No events available
    if (nEvent == 0)
        disp('ImaGIN> Warning: No events available in the file');
        return;
    % Remove all events
    elseif strcmp(EventRemove, '*')
        Event(:) = [];
        strLog = 'Removed all events';
    % Remove selected events
    else
        % Find the events corresponding to the input
        iEvt = find(strcmpi({Event.type}, EventRemove));
        % Event not found
        if isempty(iEvt)
            disp(['ImaGIN> Warning: Event "' EventRemove '" not found.']);
            return;
        end
        % Remove event
        Event(iEvt) = [];
        strLog = ['Removed event: ' EventRemove];
    end
    % Save modified events list in the file
    D = events(D,1,Event);
    save(D);
    % Add log entry
    ImaGIN_save_log(fullfile(D), strLog);
end

