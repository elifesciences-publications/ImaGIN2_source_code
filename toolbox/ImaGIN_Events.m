function ImaGIN_Events(S)

try
    t=S.Filename;
catch
    t = spm_select(1, '\.mat$', 'Select data file');
end

% Fi  = spm_figure('GetWin','Interactive');

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
  	try
        tmp=S.EventFileName{i1};
    catch
        tmp = spm_select(1, '\.txt$', sprintf('Select txt file with the timing (sec) of event %d',i1));
    end
    Timing{i1}=load(tmp);
end

try
    FileOut=S.FileOut;
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
        evt(n).time = Timing{i0}(i1);
        evt(n).value= NeventOld+i0;
    end
end
% This assigns these events to the first trials (the only one if you have
% continuous data)
D = events(D, 1, evt);
% save([pathOut 'D.filename'],D);
% P = spm_str_manip(Filename, 'H');

% 
% if str2num(version('-release'))>=14
%     save(fullfile(P, D.fname), '-V6', 'D');
% else
%     save(fullfile(P, D.fname), 'D');
% end
% 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_EventsRemove(Filename,S)

D=spm_eeg_load(Filename);

try
    OutRep=S.Out;
catch
    OutRep=[];
end

try
    EventRemove=S.EventName;
catch
    EventRemove=spm_input('Type of events to remove', '+1', 's');
end

Event=D.events;
nEvent=size(Event,2);
n=nEvent;
i1=1;
while n>0
    if strcmp(EventRemove,Event(i1).type) || isempty(EventRemove)&&isempty(Event(i1).type)
        Event(i1:nEvent-1)=Event(i1+1:nEvent);
        Event(nEvent)=[];
        i1=i1-1;
        nEvent=nEvent-1;
    end
    n=n-1;
    i1=i1+1;
end

D=events(D,1,Event);
D2=clone(D,fullfile(OutRep, D.fname), [D.nchannels D.nsamples D.ntrials]);
D2(:,:,:)=D(:,:,:);
save(D);

end


end

