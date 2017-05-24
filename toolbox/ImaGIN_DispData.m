function ImaGIN_DispData(D, SelChan)
%IMAGIN_DISPDATA Review .mat/.dat file
%
% USAGE:  ImaGIN_DispData(D,        SelChan=[Ask])
%         ImaGIN_DispData(FileName, SelChan=[Ask])
% 
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


%% ===== PARSE INPUTS =====
% Select file is not specified in input
if (nargin < 1) || isempty(D) || ~isstruct(D)
    if (nargin >= 1) && ischar(D) && ~isempty(D)
        t = D;
    else
        t = spm_select(1, '\.mat$', 'Select data file');
        if isempty(t)
            return;
        end
    end
    D = spm_eeg_load(t);
end
if (nargin < 2) || isempty(SelChan)
    SelChan = [];
end

Nchannels = nchannels(D);
Nsamples = nsamples(D);
Time = time(D);
Labels = chanlabels(D);
if (nchannels(D) == 1)
    Labels = {Labels};
end
Events = events(D);
if(~isempty(Events))
    try
        Events(1).type;
    catch
        Events = Events{1};
    end
end
Types={};
Code=[];
for i1=1:length(Events)
    trouve=0;
    for i2=1:length(Types)
        if ~strcmp(Types{i2},Events(i1).type) && i2==length(Types) && ~trouve
            Types{end+1}=Events(i1).type;
            Code(end+1)=length(Types);
        elseif strcmp(Types{i2},Events(i1).type)
            Code(end+1)=i2;
            trouve=1;
        end
    end
    if isempty(Types)
        Types{1}=Events(1).type;
        Code(1)=1;
    end
end
Ntypes=length(Types);

Data=D(:,:,:);

if isfield(D,'Firing')
    for i1=1:length(D.Firing.Name)
        Labels{end+1}=D.Firing.Name{i1};
    end
    Data=[D(:,:);D.Firing.data];
    Nchannels=Nchannels+length(D.Firing.Name);
end

if isfield(D,'Seizure')
    data=zeros(length(D.Seizure),nsamples(D));
    for i1=1:length(D.Seizure)
        data(i1,:)=D.Seizure{i1}.data;
        Labels{end+1}=D.Seizure{i1}.Name;
    end
    Data=[D(:,:,:);data];
    Nchannels=Nchannels+length(D.Seizure);
end


%% ===== CREATE FIGURE =====
FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 9);

F  = spm_figure('GetWin','Interactive');
figure(F);clf
colormap('gray')
if isempty(SelChan)
    SelChan = spm_input('Select channels', 1, 'i', sprintf('1:%d',Nchannels));
end

F  = spm_figure('GetWin','Graphics');
figure(F);clf
handles = guihandles(gcf);

M1=inf;
M2=-inf;
for i2=1:length(SelChan)
    M1=min([M1 Data(SelChan(i2),:)]);
    M2=max([M2 Data(SelChan(i2),:)]);
end
MaxNaxes=8;
% if length(SelChan)<=MaxNaxes
    for i1=1:length(SelChan)

        if mod(i1,MaxNaxes)~=0
            h.Heegaxes(i1) = axes('Position',...
                [0.13 0.83-(mod(i1,MaxNaxes)-1)*0.094 0.83 0.094],...
                'Parent', F);
        else
            h.Heegaxes(i1) = axes('Position',...
                [0.13 0.83-(MaxNaxes-1)*0.094 0.83 0.094],...
                'Parent', F);
        end
        axes(h.Heegaxes(i1));
        h.Heegplot(i1)=plot(Time,Data(SelChan(i1),:),'k');
        ylabel([sprintf('%d: ',i1) Labels{SelChan(i1)}], 'Interpreter', 'tex', 'FontSize', FS2);
        try
            axis tight
        end
        if i1==length(SelChan)|i1==MaxNaxes
            xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
        else
            set(gca,'XTickLabel',[])
        end
        set(h.Heegaxes(i1),'Visible','off')
        set(h.Heegplot(i1),'Visible','off')
        for i3=1:length(Events)
            h.eventaxes{i1,Code(i3),length(find(Code(1:i3)==Code(i3)))}=line([Events(i3).time Events(i3).time],[M1 M2],'color','m','linewidth',1,'Visible','off');
        end
    end
% end

h.DispChan=1:min([MaxNaxes length(SelChan)]);
h.c=0;
guidata(F, h);

%Display data
draw_data;





%Zoom in time
h.TimeMin=num2str(min(D.time),8);
h.TimeMax=num2str(max(D.time),8);
h.TimeCen=num2str(mean([str2double(h.TimeMin) str2double(h.TimeMax)]),4);
guidata(F, h);
nch=min([8 length(SelChan)]);
ch=nch*0.094;
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Min',...
    'Position',[0.20 0.83-ch 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Median',...
    'Position',[0.30 0.83-ch 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Max',...
    'Position',[0.40 0.83-ch 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Time',...
    'Position',[0.07 0.804-ch 0.13 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol('Tag', 'zoomtimemaxbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.40 0.804-ch 0.1 0.026],...
    'String', h.TimeMax, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomtimeminbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.20 0.804-ch 0.1 0.026],...
    'String', h.TimeMin, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomtimecenbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.30 0.804-ch 0.1 0.026],...
    'String', h.TimeCen, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomtimeinbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.19 0.052 0.035 0.026],...
    'String', '+', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,'Visible','off',...
    'Parent', F);
uicontrol('Tag', 'zoomtimeoutbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.225 0.052 0.035 0.026],...
    'String', '-', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,'Visible','off',...
    'Parent', F);
% uicontrol(F, 'Tag', 'timeslider2', 'Style', 'slider',...
%     'Min', 1, 'Max', D.nsamples, 'Value', 1, 'Units',...
%     'normalized', 'Position', [0.265 0.052 0.13 0.026],...
%     'SliderStep', [0.005 0.05],'Visible','off',...
%     'Callback', @timeslider_update,...
%     'Parent', F, 'Interruptible', 'off');

%Zoom in amplitude
h.AmpMin=num2str(M1);
h.AmpMax=num2str(M2);
h.AmpCen=num2str(mean([str2double(h.AmpMin) str2double(h.AmpMax)]));
guidata(F, h);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Amplitude',...
    'Position',[0.07 0.778-ch 0.13 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol('Tag', 'zoomampmaxbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.40 0.778-ch 0.1 0.026],...
    'String', h.AmpMax, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomampbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomampminbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.20 0.778-ch 0.1 0.026],...
    'String', h.AmpMin, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomampbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomampcenbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.30 0.778-ch 0.1 0.026],...
    'String', h.AmpCen, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomampbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomampinbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.003 0.115 0.035 0.028],...
    'String', '+', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomampbutton_update,'Visible','off',...
    'Parent', F);
uicontrol('Tag', 'zoomampoutbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.003 0.085 0.035 0.028],...
    'String', '-', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomampbutton_update,'Visible','off',...
    'Parent', F);
% uicontrol(F, 'Tag', 'freqslider2', 'Style', 'slider',...
%     'Min', 1, 'Max', M2, 'Value', 1, 'Units',...
%     'normalized', 'Position', [0.003 0.052 0.13 0.026],...
%     'SliderStep', [0.005 0.05],'Visible','off',...
%     'Callback', @ampslider_update,...
%     'Parent', F, 'Interruptible', 'off');


%Zoom
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Zoom','Value',1,'Tag','ZoomText',...
    'Position',[0.06 0.97 0.1 0.029],...
    'HorizontalAlignment', 'left', 'FontSize', FS1,...
    'BackgroundColor', 'w');
uicontrol('Tag', 'zoominbutton', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.05 0.945 0.04 0.03],...
    'String', '+', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoombutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomoutbutton', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.092 0.945 0.04 0.03],...
    'String', '-', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoombutton_update,'Visible','off',...
    'Parent', F);

%Channel slider
handles.chanrightbutton=uicontrol('Tag', 'chanrightbutton', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.945 0.05  0.04 0.03],...
    'String', '>', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @chanbutton_update,'Visible','off',...
    'Parent', F);
uicontrol('Tag', 'chanleftbutton', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.9 0.05 0.04 0.03],...
    'String', '<', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @chanbutton_update,'Visible','off',...
    'Parent', F);
if max(h.DispChan)<length(SelChan)
    set(handles.chanrightbutton,'Visible','on')
end



%Event types
if Ntypes==1
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Event type 1',...
        'Position',[0.18 0.97 0.16 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Display events',...
        'Position',[0.34 0.97 0.18 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
%     uicontrol('Tag', 'dispeventbutton', 'Style', 'radiobutton',...
%         'Units', 'normalized', 'Position', [0.41 0.945 0.04 0.03],...
%         'String', 'Y', 'FontSize', FS1,'Value',0,...
%         'BackgroundColor', spm('Colour'),...
%         'CallBack', @eventbutton_update,...
%         'Parent', F);
    uicontrol('Tag', 'dispeventbutton', 'Style', 'radiobutton',...
        'Units', 'normalized', 'Position', [0.41 0.945 0.04 0.03],...
        'Value',0,'BackgroundColor', 'w',...
        'CallBack', @eventbutton_update,...
        'Parent', F);

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Go to event',...
        'Position',[0.52 0.97 0.18 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    Neventcode=max([length(find(Code==1)) 2]);
    NeventcodeStep=max([2 Neventcode]);
    handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
        'Min', 1, 'Max', Neventcode, 'Value', 1, 'Units',...
        'normalized', 'Position', [0.53 0.957 0.16 0.02],...
        'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 10/(NeventcodeStep-1))],...
        'TooltipString', 'Choose trial number',...
        'Callback', @eventslider_update,...
        'Parent', F, 'Interruptible', 'off');
    % frame for eventgoslider text
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.53 0.932 0.16 0.025]);
    % event go slider texts
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', '1',...
        'Position',[0.535 0.933 0.03 0.022],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    handles.eventgotext = uicontrol(F, 'Style', 'text', 'Tag', 'eventgotext',...
        'Units', 'normalized',...
        'String', int2str(get(handles.eventgoslider, 'Value')),...
        'Position',[0.59 0.933 0.04 0.022],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    uicontrol(F,'Tag', 'eventgotextmax', 'Style', 'text', 'Units', 'normalized',...
        'String', num2str(Neventcode),...
        'Position',[0.645 0.933 0.04 0.022],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Event type name',...
        'Position',[0.74 0.97 0.20 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.74 0.95 0.20 0.025]);
    uicontrol(F, 'Style', 'text', ...
        'Units', 'normalized',...
        'String', Events(1).type,...
        'Position',[0.745 0.951 0.19 0.022],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));

    handles.eventtypetext = uicontrol(F, 'Style', 'text', 'Tag', 'eventtypetext',...
        'Units', 'normalized',...
        'String', '1',...
        'Position',[0.24 0.933 0.04 0.022],'Visible','off',...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));


elseif Ntypes>1
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Event type',...
        'Position',[0.18 0.97 0.16 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    handles.eventtypeslider = uicontrol(F, 'Tag', 'eventtypeslider', 'Style', 'slider',...
        'Min', 1, 'Max', Ntypes, 'Value', 1, 'Units',...
        'normalized', 'Position', [0.18 0.957 0.16 0.02],...
        'SliderStep', [1/(Ntypes-1) min(Ntypes-1, 10/(Ntypes-1))],...
        'TooltipString', 'Choose trial number',...
        'Callback', @eventtypeslider_update,...
        'Parent', F, 'Interruptible', 'off');
    % frame for eventtypeslider text
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.18 0.932 0.16 0.025]);
    % event types slider texts
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', '1',...
        'Position',[0.185 0.933 0.03 0.022],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Display events',...
        'Position',[0.34 0.97 0.18 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    handles.eventtypetext = uicontrol(F, 'Style', 'text', 'Tag', 'eventtypetext',...
        'Units', 'normalized',...
        'String', int2str(get(handles.eventtypeslider, 'Value')),...
        'Position',[0.24 0.933 0.04 0.022],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', mat2str(Ntypes),...
        'Position',[0.29 0.933 0.04 0.022],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    uicontrol('Tag', 'dispeventbutton', 'Style', 'radiobutton',...
        'Units', 'normalized', 'Position', [0.41 0.945 0.04 0.03],...
        'FontSize', FS1,'Value',0,...
        'BackgroundColor', 'w',...
        'CallBack', @eventbutton_update,...
        'Parent', F);

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Go to event',...
        'Position',[0.52 0.97 0.18 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    Neventcode=max([length(find(Code==1)) 2]);
    NeventcodeStep=max([2 Neventcode]);
    handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
        'Min', 1, 'Max', Neventcode, 'Value', 1, 'Units',...
        'normalized', 'Position', [0.53 0.957 0.16 0.02],...
        'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 10/(NeventcodeStep-1))],...
        'TooltipString', 'Choose trial number',...
        'Callback', @eventslider_update,...
        'Parent', F, 'Interruptible', 'off');
    % frame for eventgoslider text
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.53 0.932 0.16 0.025]);
    % event go slider texts
    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', '1',...
        'Position',[0.535 0.933 0.03 0.022],...
        'HorizontalAlignment', 'left', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    handles.eventgotext = uicontrol(F, 'Style', 'text', 'Tag', 'eventgotext',...
        'Units', 'normalized',...
        'String', int2str(get(handles.eventgoslider, 'Value')),...
        'Position',[0.59 0.933 0.04 0.022],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    uicontrol(F,'Tag', 'eventgotextmax', 'Style', 'text', 'Units', 'normalized',...
        'String', num2str(Neventcode),...
        'Position',[0.645 0.933 0.04 0.022],...
        'HorizontalAlignment', 'right', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));

    uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Event type name',...
        'Position',[0.74 0.97 0.20 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
        'normalized', 'Position',[0.74 0.95 0.20 0.025]);
    uicontrol(F, 'Style', 'text','Tag','eventnametext',...
        'Units', 'normalized',...
        'String', Types{get(handles.eventtypeslider, 'Value')},...
        'Position',[0.745 0.951 0.19 0.022],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
end

handles.timeslider=uicontrol(F, 'Tag', 'timeslider', 'Style', 'slider',...
    'Min', 1, 'Max', Nsamples, 'Value', 1, 'Units',...
    'normalized', 'Position', [0.625 0.82-ch 0.2 0.03],...
    'SliderStep', [0.005 0.05],...
    'TooltipString', 'Choose sample number',...
    'Callback', @timeslider_update,...
    'Parent', F, 'Interruptible', 'off');
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Time',...
    'Position',[0.625 0.776-ch 0.2 0.029],...
    'HorizontalAlignment', 'center', 'FontSize', FS1,...
    'BackgroundColor', 'w');


h.D=D;
h.Data=Data;
h.Events=Events;
h.SelChan=SelChan;
h.MaxNaxes=MaxNaxes;
h.M1=M1;
h.M2=M2;
h.offset=1;
h.Time=Time;
h.Types=Types;
h.Ntypes=Ntypes;
h.Code=Code;
h.ZoomAmp=0;
h.zoomb=0;
% store handles
guidata(F, h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomampbutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
value = str2num(get(hObject,'String'));
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
M1=inf;
M2=-inf;
for i2=1:length(h.SelChan)
    M1=min([M1 h.Data(h.SelChan(i2),:)]);
    M2=max([M2 h.Data(h.SelChan(i2),:)]);
end
switch action
    case 'zoomampmaxbutton'
        h.AmpMax=min([value M2]);
        h.AmpMax=num2str(h.AmpMax);
        h.AmpCen=num2str(mean([str2num(h.AmpMin) str2num(h.AmpMax)]));
        set(handles.zoomampcenbutton,'String',h.AmpCen)
        set(handles.zoomampmaxbutton,'String',h.AmpMax)
    case 'zoomampminbutton'
        h.AmpMin=max([value M1]);
        h.AmpMin=num2str(h.AmpMin);
        h.AmpCen=num2str(mean([str2num(h.AmpMin) str2num(h.AmpMax)]));
        set(handles.zoomampminbutton,'String',h.AmpMin)
        set(handles.zoomampcenbutton,'String',h.AmpCen)
    case 'zoomampcenbutton'
        h.AmpMin=num2str(str2num(h.AmpMin)+(value-str2num(h.AmpCen)));
        set(handles.zoomampminbutton,'String',h.AmpMin)
        h.AmpMax=num2str(str2num(h.AmpMax)+(value-str2num(h.AmpCen)));
        set(handles.zoomampmaxbutton,'String',h.AmpMax)
        h.AmpCen=num2str(value);
end

z=get(handles.ZoomText,'Value');
offset=h.offset;
h.ZoomAmp=1;
guidata(F,h)


% make plots
if strcmp(action,'zoomampmaxbutton')|strcmp(action,'zoomampminbutton')|strcmp(action,'zoomampcenbutton')
     draw_subplots(handles, z, offset)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomtimebutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
value = str2num(get(hObject,'String'));
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
D=h.D;

switch action
    case 'zoomtimemaxbutton'
        h.TimeMax=min([D.time(end) value]);
        h.TimeMax=num2str(h.TimeMax);
        h.TimeCen=num2str(mean([str2num(h.TimeMin) str2num(h.TimeMax)]));
        set(handles.zoomtimecenbutton,'String',h.TimeCen)
        set(handles.zoomtimemaxbutton,'String',h.TimeMax)
    case 'zoomtimeminbutton'
        h.TimeMin=max([D.time(1) value]);
        h.TimeMin=num2str(h.TimeMin);
        h.TimeCen=num2str(mean([str2num(h.TimeMin) str2num(h.TimeMax)]));
        set(handles.zoomtimecenbutton,'String',h.TimeCen)
        set(handles.zoomtimeminbutton,'String',h.TimeMin)
    case 'zoomtimecenbutton'
        h.TimeMin=num2str(str2num(h.TimeMin)+(value-str2num(h.TimeCen)));
        set(handles.zoomtimeminbutton,'String',h.TimeMin)
        h.TimeMax=num2str(str2num(h.TimeMax)+(value-str2num(h.TimeCen)));
        set(handles.zoomtimemaxbutton,'String',h.TimeMax)
        h.TimeCen=num2str(value);
end
guidata(F,h)
samplemin=indsample(D,str2num(h.TimeMin));
samplemax=indsample(D,str2num(h.TimeMax));
%update time slider
z=D.nsamples/(samplemax-samplemin);
offset= samplemin+round(0.4*D.nsamples/z);
set(handles.timeslider, 'Value',offset);
set(handles.ZoomText, 'Value', z);
h.offset=offset;


% make plots
if strcmp(action,'zoomtimemaxbutton')|strcmp(action,'zoomtimeminbutton')|strcmp(action,'zoomtimecenbutton')
    draw_subplots(handles, z, offset)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoombutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

offset=h.offset;

if strcmp(action,'zoominbutton')
    z=get(handles.ZoomText,'Value');
    z=z*1.2;
    set(handles.ZoomText,'Value',z)
elseif strcmp(action,'zoomoutbutton')
    z=get(handles.ZoomText,'Value');
    if z>1
        z=z/1.2;
        set(handles.ZoomText,'Value',z)
    end
end
if z>1
    set(handles.zoomoutbutton,'Visible','on')
else
    set(handles.zoomoutbutton,'Visible','off')
end

if z<1
    z=1;
    set(handles.ZoomText,'Value',z)
end
h.zoomb=1;
guidata(F,h);

% make plots
draw_subplots(handles, z,offset)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chanbutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

% if strcmp(action,'chanrightbutton')
%     h.c=h.c+1;
%     h.DispChan=h.c*h.MaxNaxes+[1:h.MaxNaxes];
%     h.DispChan=h.DispChan(find(h.DispChan<=length(h.SelChan)));
% %     h.DispChan=max(h.DispChan)+[1:h.MaxNaxes];
%     if max(h.DispChan)>=length(h.SelChan)
%         set(handles.chanrightbutton,'Visible','off')
% %         h.DispChan=[max(h.SelChan)-h.MaxNaxes+1:max(h.SelChan)];
%     end
%     if min(h.DispChan)>1
%         set(handles.chanleftbutton,'Visible','on')
%     else
%         set(handles.chanleftbutton,'Visible','off')
%     end
% elseif strcmp(action,'chanleftbutton')
%     h.c=h.c-1;
%     h.DispChan=[min(h.DispChan)-h.MaxNaxes:min(h.DispChan)-1];
%     if min(h.DispChan)<=1
%         set(handles.chanleftbutton,'Visible','off')
%         h.DispChan=[1:h.MaxNaxes];
%     end
%     if max(h.DispChan)<length(h.SelChan)
%         set(handles.chanrightbutton,'Visible','on')
%     else
%         set(handles.chanrightbutton,'Visible','off')
%     end
% end

if strcmp(action,'chanrightbutton')
    h.c=h.c+1;
elseif strcmp(action,'chanleftbutton')
    h.c=h.c-1;
end
h.DispChan=h.c*h.MaxNaxes+[1:h.MaxNaxes];
h.DispChan=h.DispChan(find(h.DispChan<=length(h.SelChan)&h.DispChan>=1));
if max(h.DispChan)>=length(h.SelChan)
    set(handles.chanrightbutton,'Visible','off')
else
    set(handles.chanrightbutton,'Visible','on')
end
if min(h.DispChan)>1
    set(handles.chanleftbutton,'Visible','on')
else
    set(handles.chanleftbutton,'Visible','off')
end

guidata(F,h);

% make plots
draw_data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventbutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Value');

% if action
%     set(hObject,'Value',0,'String','Y');
% else
%     set(hObject,'Value',1,'String','N');
% end
    
% make plots
draw_bars(action)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventtypeslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
offset=h.offset;
Events=h.Events;
Code=h.Code;

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

% Update display of current trial number
set(handles.eventtypetext, 'String', mat2str(ind));
set(handles.eventnametext, 'String', h.Types{ind})

% Update slider for events
Neventcode=max([length(find(Code==ind)) 2]);
NeventcodeStep=max([2 Neventcode]);
set(handles.eventgoslider,'Max', Neventcode, 'Value', 1,'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 10/(NeventcodeStep-1))]);
set(handles.eventgotext, 'String', '1');
set(handles.eventgotextmax, 'String', num2str(Neventcode))



% make plots
draw_subplots(handles, get(handles.ZoomText,'Value'),offset)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
Code=h.Code;

ntype=str2num(get(handles.eventtypetext, 'String'));

% slider value
ind = round(get(hObject, 'Value'));
tmp=length(find(Code==ntype));
if ind>tmp
    ind=tmp;
end
set(hObject, 'Value', ind);

% Update display of current trial number
set(handles.eventgotext, 'String', mat2str(ind));

nevent=str2num(get(handles.eventgotext, 'String'));
tmp=find(Code==ntype);
Events=h.Events;
offset=indsample(D,Events(tmp(nevent)).time);
set(handles.timeslider, 'Value', offset);
h.offset=offset;
guidata(F,h);


% make plots
draw_subplots(handles, get(handles.ZoomText,'Value'), offset)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

offset=ind;
z=get(handles.ZoomText,'Value');

h.offset=offset;
h.zoomb=1;
guidata(F,h);

% make plots
draw_subplots(handles, z, offset)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_subplots(handles, z, offset)
% This function plots data in the subplots

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

D=h.D;
Data=h.Data;
Events=h.Events;
Code=h.Code;
SelChan=h.SelChan;
MaxNaxes=h.MaxNaxes;
M1=h.M1;
M2=h.M2;

if h.Ntypes>0
    Ntype=str2num(get(handles.eventtypetext,'String'));
    Ngo=str2num(get(handles.eventgotext,'String'));
    n=length(find(h.Code==Ntype));
else
    n=0;
end

Time=h.Time;


TimeIndex=offset+[1:round(D.nsamples/z)]-round(0.4*round(D.nsamples/z));

if min(TimeIndex)<1
    TimeIndex=1-min(TimeIndex)+TimeIndex;
elseif max(TimeIndex)>D.nsamples
    TimeIndex=D.nsamples-max(TimeIndex)+TimeIndex;
end

if h.zoomb
    %update zoomampbutton
    h.TimeMin=num2str(min(TimeIndex/D.fsample));
    h.TimeMax=num2str(max(TimeIndex/D.fsample));
    h.TimeCen=num2str(mean([h.TimeMin h.TimeMax]));
    set(handles.zoomtimeminbutton,'String',h.TimeMin);
    set(handles.zoomtimemaxbutton,'String',h.TimeMax);
    set(handles.zoomtimecenbutton,'String',h.TimeCen);
end

for i1=1:length(SelChan)

        axes(h.Heegaxes(i1))
       % try
            if h.ZoomAmp
                axis([min(Time(TimeIndex)) max(Time(TimeIndex)) str2num(h.AmpMin) str2num(h.AmpMax)])
            else
                axis([min(Time(TimeIndex)) max(Time(TimeIndex)) min(Data(SelChan(i1),TimeIndex)) max(Data(SelChan(i1),TimeIndex))])
            end
       % end
        set(h.Heegaxes(i1),'Visible','off')
        set(h.Heegplot(i1),'Visible','off')

        for i3=1:length(Events)
            set(h.eventaxes{i1,Code(i3),length(find(Code(1:i3)==Code(i3)))},'Visible','off')
            h.eventaxes{i1,Code(i3),length(find(Code(1:i3)==Code(i3)))}=line([Events(i3).time Events(i3).time],[M1 M2],'color','m','linewidth',1,'Visible','off');
        end

        for i2=1:n
            if get(handles.dispeventbutton,'Value')
                set(h.eventaxes{i1,Ntype,i2},'Visible','on')
%             else
%                 set(h.eventaxes{i1,Ntype,i2},'Visible','off')
            end
        end

    
    end
    draw_data;
    h.ZoomAmp=0;
    h.zoomb=0;
guidata(F,h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_bars(action)
% This function plots events in the subplots

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
Types=h.Types;
Nevents=length(Types);
Events=h.Events;

if Nevents>1
    Ntype=str2num(get(handles.eventtypetext,'String'));
else
    Ntype=1;
end
n=0;
for i1=1:length(Events)
    if strcmp(Events(i1).type,Types{Ntype})
        n=n+1;
        Eventref(n)=i1;
    end
end
n=length(Eventref);

SelChan=h.SelChan;
MaxNaxes=h.MaxNaxes;
% if length(SelChan)<=MaxNaxes
    for i1=1:length(SelChan)
        for i2=1:n
            if action
                set(h.eventaxes{i1,Ntype,i2},'Visible','on')
            else
                set(h.eventaxes{i1,Ntype,i2},'Visible','off')
            end
            tmp=setdiff(Nevents,Ntype);
            for i3=1:length(tmp)
                set(h.eventaxes{i1,tmp(i3),i2},'Visible','off')
            end
        end
    end
% end
guidata(F,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_data
% This function plots data in the subplots

FS1 = spm('FontSize', 14);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

for i1=1:length(h.Heegaxes)
    if isempty(intersect(i1,h.DispChan))
        set(h.Heegaxes(i1),'Visible','off')
        set(h.Heegplot(i1),'Visible','off')
    else
        set(h.Heegaxes(i1),'Visible','on')
        set(h.Heegplot(i1),'Visible','on')
        axes(h.Heegaxes(i1))
        if i1==max(h.DispChan)
            xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
            set(h.Heegaxes(i1),'XTickLabelMode','auto')
        else
            xlabel('');
%             set(h.Heegaxes(i1),'XTickLabel',[])
        end
    end
end
guidata(F,h)

