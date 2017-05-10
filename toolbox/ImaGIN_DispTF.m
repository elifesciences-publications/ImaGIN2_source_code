function ImaGIN_DispTF(varargin)
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

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 10);

t = spm_select(1, '\.mat$', 'Select data file');
D = spm_eeg_load(t);

Events = events(D);
if ~isempty(Events)
    try
        Events(1).type;
    catch
        Events=Events{1};
    end
end
Types = {};
Code = [];
for i1 = 1:length(Events)
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
Ntypes = length(Types);

% Load and clear interactive window
Fi  = spm_figure('GetWin','Interactive');
figure(Fi);
clf;

if ~isfield(D,'tf')
    Ctype = {
        'Morlet wavelet',...
        'Hilbert',...
        'Coherence'};
    str   = 'Time frequency decomposition';
    Sel   = spm_input(str, 2, 'm', Ctype);
    Method = Ctype{Sel};

    Flag = spm_input('Display ',1,'Power|Synchrony');

    Direc = spm_str_manip(t,'h');
    File = spm_str_manip(t,'t');
    E = what(Direc);
    ok = 0;
    switch Method
        case 'Morlet wavelet'
            % Check if Morlet wavelet has been computed
            for i1=1:length(E.mat)
                if strcmp(E.mat{i1},['w1_' File])
                    ok=1;
                    break
                end
            end
            if ~ok
                spm_input('Morlet wavelet not computed yet', '+1','d')
                return
            else
                switch Flag
                    case{'Power'}
                        Flag2 = spm_input('Power integrated over time window ','+1','Yes|No');
                        if strcmp(Flag2,'Yes')
                            t=fullfile(Direc,['w1_' File]);
                        else
                            t=fullfile(Direc,['w1_' File]);
                        end
                    case{'Synchrony'}
                        t=fullfile(Direc,['w2_' File]);
                end
            end

        case 'Hilbert'
            % Check if Hilbert has been computed
            for i1=1:length(E.mat)
                if strcmp(E.mat{i1},['h1_' File])
                    ok=1;
                    break
                end
            end
            if ~ok
                spm_input('Hilbert not computed yet', '+1','d')
                return
            else
                switch Flag
                    case 'Power'
                        Flag2 = spm_input('Power integrated over time window ','+1','Yes|No');
                        if strcmp(Flag2,'Yes')
                            t=fullfile(Direc,['h1_' File]);
                        else
                            t=fullfile(Direc,['h1_' File]);
                        end
                    case 'Synchrony'
                        t=fullfile(Direc,['h2_' File]);
                end
            end

        case 'Coherence'
            % Check if coherence has been computed
            for i1=1:length(E.mat)
                if strcmp(E.mat{i1},['c1_' File])
                    ok = 1;
                    break;
                end
            end
            if ~ok
                spm_input('Coherence not computed yet', '+1','d')
                return;
            else
                switch Flag
                    case 'Power'
                        t = fullfile(Direc,['c1_' File]);
                    case 'Synchrony'
                        t = fullfile(Direc,['c2_' File]);
                end
            end
    end
    % Load data
    D = spm_eeg_load(t);

else
    Method = D.tf.Method;
end
Flag = D.tf.Label;

% Reinterpolate for display
Freq = min(D.tf.frequencies):min([(max(D.tf.frequencies)-min(D.tf.frequencies))/300 min(diff(D.tf.frequencies))]):max(D.tf.frequencies);
Temps = min(D.tf.time):min([(max(D.tf.time)-min(D.tf.time))/1e3 min(diff(D.tf.time))]):max(D.tf.time);
ok = 1;
while ok
    try
        ok=0;
        DataDisp=zeros(D.nchannels,length(Freq),length(Temps));
    catch
        ok=1;
        Freq=Freq(1:2:end);
        Temps=Temps(1:2:end);
    end
end
    
[x,y]=meshgrid(D.tf.time(1:size(D,3)),D.tf.frequencies);
[xi,yi]=meshgrid(Temps,Freq);

for i1=1:D.nchannels
    for k=1:D.ntrials
        DataDisp(i1,:,:,k)=interp2(x,y,squeeze(D(i1,:,:,k)),xi,yi);
    end
end
D.tf.time=Temps;
D.tf.frequencies=Freq;

% Load Graphics window
F  = spm_figure('GetWin','Graphics');
figure(F);clf
handles = guihandles(gcf);

handles.dispeventbutton=uicontrol('Tag', 'dispeventbutton', 'Style', 'radiobutton',...
    'Units', 'normalized', 'Position', [0.41 0.945 0.04 0.03],...
    'FontSize', FS1,'Value',0,'Visible','off',...
    'BackgroundColor', 'w',...
    'CallBack', @eventbutton_update,...
    'Parent', F);

handles.normalisebutton=uicontrol('Tag', 'normalisebutton', 'Style', 'radiobutton',...
    'Units', 'normalized', 'Position', [0.95 0.484 0.04 0.03],...
    'FontSize', FS1,'Value',0,'Visible','off',...
    'BackgroundColor', 'w',...
    'CallBack', @normalisebutton_update,...
    'Parent', F);


% All channel selected by default
SelChan=1:D.nchannels;
SelTrial=1:D.ntrials;

% Store data
h.Flag = Flag;
try
    time=D.tf.time;
catch
    time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
    time=time-time(D.TimeZero);
    D.tf.time=time;
end
h.time=time;
h.D=D;
h.Events=Events;
h.ntype=1;
h.SelChan=SelChan;
h.SelTrial=SelTrial;
h.offsettime=1;
h.offsetfreq=1;
h.ztime=1;
h.zfreq=1;
h.zcol=1;
h.offsettime2=1;
h.offsetfreq2=1;
h.ztime2=1;
h.zfreq2=1;
h.zcol2=1;
h.Method=Method;
h.Types=Types;
h.Ntypes=Ntypes;
h.Code=Code;
h.DataDisp=DataDisp;
h.poweraxes = axes('Position',...
    [0.09 0.55 0.87 0.35],...
    'Parent', F);
guidata(F, h);

% Compute power averaged on electrodes
data_plots(h.DataDisp,1);

% Make plots without bars
draw_plots(1,0)

h = guidata(F);

% Zoom in time
h.TimeMin=num2str(min(D.tf.time),4);
h.TimeMax=num2str(max(D.tf.time),4);
h.TimeCen=num2str(mean([str2double(h.TimeMin) str2double(h.TimeMax)]),4);
guidata(F, h);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Min',...
    'Position',[0.14 0.481 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Median',...
    'Position',[0.24 0.481 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Max',...
    'Position',[0.34 0.481 0.1 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Time',...
    'Position',[0.01 0.455 0.13 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol('Tag', 'zoomtimemaxbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.34 0.455 0.1 0.026],...
    'String', h.TimeMax, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomtimeminbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.14 0.455 0.1 0.026],...
    'String', h.TimeMin, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomtimebutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomtimecenbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.24 0.455 0.1 0.026],...
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
uicontrol(F, 'Tag', 'timeslider2', 'Style', 'slider',...
    'Min', 1, 'Max', D.nsamples, 'Value', 1, 'Units',...
    'normalized', 'Position', [0.265 0.052 0.13 0.026],...
    'SliderStep', [0.005 0.05],'Visible','off',...
    'Callback', @timeslider_update,...
    'Parent', F, 'Interruptible', 'off');

% Zoom in frequency
h.FreqMin=num2str(min(D.tf.frequencies));
h.FreqMax=num2str(max(D.tf.frequencies));
h.FreqCen=num2str(mean([str2double(h.FreqMin) str2double(h.FreqMax)]));
guidata(F, h);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Frequency',...
    'Position',[0.01 0.429 0.13 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol('Tag', 'zoomfreqmaxbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.34 0.429 0.1 0.026],...
    'String', h.FreqMax, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomfreqbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomfreqminbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.14 0.429 0.1 0.026],...
    'String', h.FreqMin, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomfreqbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomfreqcenbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.24 0.429 0.1 0.026],...
    'String', h.FreqCen, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomfreqbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomfreqinbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.003 0.115 0.035 0.028],...
    'String', '+', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomfreqbutton_update,'Visible','off',...
    'Parent', F);
uicontrol('Tag', 'zoomfreqoutbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.003 0.085 0.035 0.028],...
    'String', '-', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomfreqbutton_update,'Visible','off',...
    'Parent', F);
uicontrol(F, 'Tag', 'freqslider2', 'Style', 'slider',...
    'Min', 1, 'Max', D.Nfrequencies, 'Value', 1, 'Units',...
    'normalized', 'Position', [0.003 0.052 0.13 0.026],...
    'SliderStep', [0.005 0.05],'Visible','off',...
    'Callback', @freqslider_update,...
    'Parent', F, 'Interruptible', 'off');

% Frequency normalisation
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Frequency normalisation',...
    'Position',[0.65 0.482 0.3 0.029],...
    'HorizontalAlignment', 'center', 'FontSize', FS1,...
    'BackgroundColor', 'w',...
    'Visible','off');
set(handles.normalisebutton,'Visible','off');


% Adjust colours
h.ColMin=num2str(h.MinPlot,4);
h.ColMax=num2str(h.MaxPlot,4);
h.ColCen=num2str(mean([str2num(h.ColMin) str2num(h.ColMax)]),4);
guidata(F, h);
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Colour',...
    'Position',[0.01 0.403 0.13 0.026],...
    'BackgroundColor', spm('Colour'),...
    'HorizontalAlignment', 'center', 'FontSize', FS2);
uicontrol('Tag', 'zoomcolmaxbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.34 0.403 0.1 0.026],...
    'String', h.ColMax, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomcolbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomcolminbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.14 0.403 0.1 0.026],...
    'String', h.ColMin, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomcolbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomcolcenbutton', 'Style', 'edit',...
    'Units', 'normalized', 'Position', [0.24 0.403 0.1 0.026],...
    'String', h.ColCen, 'FontSize', FS3,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomcolbutton_update,...
    'Parent', F);
uicontrol('Tag', 'zoomcolinbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.965 0.275 0.035 0.028],...
    'String', '+', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomcolbutton_update,'Visible','off',...
    'Parent', F);
uicontrol('Tag', 'zoomcoloutbutton2', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.965 0.245 0.035 0.028],...
    'String', '-', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @zoomcolbutton_update,'Visible','off',...
    'Parent', F);

%Event-related averaging
uicontrol('Tag', 'averagebutton', 'Style', 'pushbutton',...
    'Units', 'normalized', 'Position', [0.725 0.915 0.23 0.029],...
    'String', 'Event-related average', 'FontSize', FS2,...
    'BackgroundColor', spm('Colour'),...
    'CallBack', @averagebutton,...
    'Parent', F,'visible','off');



%Channel Selection
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'Channel',...
    'Position',[0.01 0.97 0.16 0.029],...
    'HorizontalAlignment', 'center', 'FontSize', FS1,...
    'BackgroundColor', 'w');
handles.channelslider = uicontrol(F, 'Tag', 'channelslider', 'Style', 'slider',...
    'Min', 1, 'Max', D.nchannels+1, 'Value', 1, 'Units',...
    'normalized', 'Position', [0.01 0.957 0.16 0.02],...
    'SliderStep', [1/(D.nchannels) min(D.nchannels, 1/(D.nchannels))],...
    'Callback', @channelslider_update,...
    'Parent', F, 'Interruptible', 'off');
% frame for channelslider text
uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
    'normalized', 'Position',[0.01 0.932 0.16 0.025]);
% channel slider texts
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', 'All',...
    'Position',[0.012 0.933 0.03 0.022],...
    'HorizontalAlignment', 'left', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'));
uicontrol(F, 'Style', 'text', 'Tag', 'channeltext',...
    'Units', 'normalized',...
    'String', 'All',...
    'Position',[0.057 0.933 0.05 0.022],...
    'HorizontalAlignment', 'center', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'));
uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
    'String', mat2str(D.nchannels),...
    'Position',[0.115 0.933 0.05 0.022],...
    'HorizontalAlignment', 'right', 'FontSize', FS1,...
    'BackgroundColor', spm('Colour'));


%Determine if several trials
FlagTrial=0;
if size(D,4)>1
    FlagTrial=1;    %then events are trials...
end
                
if FlagTrial
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
        set(handles.dispeventbutton,'Visible','on');

        Neventcode=D.ntrials;

        if Neventcode>1
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to event',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w');
            handles.trialgoslider = uicontrol(F, 'Tag', 'trialgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', Neventcode+1, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.02],...
                'SliderStep', [1/(Neventcode) min(Neventcode, 1/(Neventcode))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @trialslider_update,...
                'Parent', F, 'Interruptible', 'off');
            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','on');
            % event go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'All',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            handles.trialgotext = uicontrol(F, 'Style', 'text', 'Tag', 'trialgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.trialgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
        else
            NeventcodeStep=2;
            handles.trialgoslider = uicontrol(F, 'Tag', 'trialgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', NeventcodeStep, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.022],...
                'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @trialslider_update,...
                'Parent', F, 'Interruptible', 'off','Visible','off');

            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','off');
            % event go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', '1',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            handles.trialgotext = uicontrol(F, 'Style', 'text', 'Tag', 'trialgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.trialgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
        end

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
            'SliderStep', [1/(Ntypes-1) min(Ntypes-1, 1/(Ntypes-1))],...
            'TooltipString', 'Choose event type',...
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
            'BackgroundColor', 'w');
        
        set(handles.dispeventbutton,'Visible','on');

        Neventcode=D.ntrials;
        if Neventcode>1
            handles.trialgoslider = uicontrol(F, 'Tag', 'trialgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', Neventcode+1, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.022],...
                'SliderStep', [1/(Neventcode) min(Neventcode, 1/(Neventcode))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @trialslider_update,...
                'Parent', F, 'Interruptible', 'off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to trial',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w');
            % frame for trialgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','on');
            % trial go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'All',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            handles.trialgotext = uicontrol(F, 'Style', 'text', 'Tag', 'trialgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.trialgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
        else
            NeventcodeStep=2;
            handles.trialgoslider = uicontrol(F, 'Tag', 'trialgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', NeventcodeStep, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.022],...
                'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @trialslider_update,...
                'Parent', F, 'Interruptible', 'off','Visible','off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to trial',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w','Visible','off');
            % frame for trialgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','off');
            % trial go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', '1',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            handles.trialgotext = uicontrol(F, 'Style', 'text', 'Tag', 'trialgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.trialgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
        end

        uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
            'String', 'Event type name',...
            'Position',[0.74 0.97 0.20 0.029],...
            'HorizontalAlignment', 'center', 'FontSize', FS1,...
            'BackgroundColor', 'w');
        uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
            'normalized', 'Position',[0.74 0.95 0.20 0.025]);
        uicontrol(F, 'Style', 'text','Tag','eventnametext',...
            'Units', 'normalized',...
            'String', Events(get(handles.eventtypeslider, 'Value')).type,...
            'Position',[0.745 0.951 0.19 0.022],...
            'HorizontalAlignment', 'center', 'FontSize', FS1,...
            'BackgroundColor', spm('Colour'));

    end
else
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

        set(handles.dispeventbutton,'Visible','on');

        Neventcode=length(find(Code==1));

        if Neventcode>1
            handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', Neventcode, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.02],...
                'SliderStep', [1/(Neventcode-1) min(Neventcode-1, 1/(Neventcode-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @eventslider_update,...
                'Parent', F, 'Interruptible', 'off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to event',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w');
            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','on');
            % event go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', '1',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            handles.eventgotext = uicontrol(F, 'Style', 'text', 'Tag', 'eventgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.eventgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
        else
            NeventcodeStep=2;
            handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', NeventcodeStep, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.02],...
                'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @eventslider_update,...
                'Parent', F, 'Interruptible', 'off','Visible','off');

            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','off');
            % event go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', '1',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            handles.eventgotext = uicontrol(F, 'Style', 'text', 'Tag', 'eventgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.eventgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','off');
        end

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
            'SliderStep', [1/(Ntypes-1) min(Ntypes-1, 1/(Ntypes-1))],...
            'TooltipString', 'Choose event number',...
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
        set(handles.dispeventbutton,'Visible','on');

         Neventcode=max([length(Events(1).time) 2]);
         NeventcodeStep=max([2 Neventcode]);
        if Neventcode>1
            handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', Neventcode, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.02],...
                'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @eventslider_update,...
                'Parent', F, 'Interruptible', 'off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to event',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w');
            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','on');
            % event go slider texts
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', '1',...
                'Position',[0.535 0.933 0.03 0.022],...
                'HorizontalAlignment', 'left', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            handles.eventgotext = uicontrol(F, 'Style', 'text', 'Tag', 'eventgotext',...
                'Units', 'normalized',...
                'String', int2str(get(handles.eventgoslider, 'Value')),...
                'Position',[0.59 0.933 0.04 0.022],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', num2str(Neventcode),...
                'Position',[0.645 0.933 0.04 0.022],...
                'HorizontalAlignment', 'right', 'FontSize', FS1,...
                'BackgroundColor', spm('Colour'),'Visible','on');
        else
            NeventcodeStep=2;
            handles.eventgoslider = uicontrol(F, 'Tag', 'eventgoslider', 'Style', 'slider',...
                'Min', 1, 'Max', NeventcodeStep, 'Value', 1, 'Units',...
                'normalized', 'Position', [0.53 0.957 0.16 0.02],...
                'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))],...
                'TooltipString', 'Choose trial number',...
                'Callback', @eventslider_update,...
                'Parent', F, 'Interruptible', 'off','Visible','off');
            uicontrol(F, 'Style', 'text', 'Units', 'normalized',...
                'String', 'Go to event',...
                'Position',[0.52 0.97 0.18 0.029],...
                'HorizontalAlignment', 'center', 'FontSize', FS1,...
                'BackgroundColor', 'w','Visible','off');
            % frame for eventgoslider text
            uicontrol(F, 'Style','Frame','BackgroundColor',spm('Colour'), 'Units',...
                'normalized', 'Position',[0.53 0.932 0.16 0.025],'Visible','off');
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
        end

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
end
guidata(F, h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomfreqbutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
value = str2num(get(hObject,'String'));
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

switch action
    case 'zoomfreqmaxbutton'
        h.FreqMax=num2str(value);
        h.FreqCen=num2str(mean([str2num(h.FreqMin) str2num(h.FreqMax)]));
        set(handles.zoomfreqcenbutton,'String',h.FreqCen)
    case 'zoomfreqminbutton'
        h.FreqMin=num2str(value);
        h.FreqCen=num2str(mean([str2num(h.FreqMin) str2num(h.FreqMax)]));
        set(handles.zoomfreqcenbutton,'String',h.FreqCen)
    case 'zoomfreqcenbutton'
        h.FreqMin=num2str(str2num(h.FreqMin)+(value-str2num(h.FreqCen)));
        set(handles.zoomfreqminbutton,'String',h.FreqMin)
        h.FreqMax=num2str(str2num(h.FreqMax)+(value-str2num(h.FreqCen)));
        set(handles.zoomfreqmaxbutton,'String',h.FreqMax)
        h.FreqCen=num2str(value);
end
guidata(F,h)

% make plots
if strcmp(action,'zoomfreqmaxbutton') || strcmp(action,'zoomfreqminbutton') || strcmp(action,'zoomfreqcenbutton')
    adjust_window(h.poweraxes)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomtimebutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
value = str2num(get(hObject,'String'));
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

switch action
    case 'zoomtimemaxbutton'
        h.TimeMax=num2str(value);
        h.TimeCen=num2str(mean([str2num(h.TimeMin) str2num(h.TimeMax)]));
        set(handles.zoomtimecenbutton,'String',h.TimeCen)
    case 'zoomtimeminbutton'
        h.TimeMin=num2str(value);
        h.TimeCen=num2str(mean([str2num(h.TimeMin) str2num(h.TimeMax)]));
        set(handles.zoomtimecenbutton,'String',h.TimeCen)
    case 'zoomtimecenbutton'
        h.TimeMin=num2str(str2num(h.TimeMin)+(value-str2num(h.TimeCen)));
        set(handles.zoomtimeminbutton,'String',h.TimeMin)
        h.TimeMax=num2str(str2num(h.TimeMax)+(value-str2num(h.TimeCen)));
        set(handles.zoomtimemaxbutton,'String',h.TimeMax)
        h.TimeCen=num2str(value);
end
guidata(F,h)

% make plots
if strcmp(action,'zoomtimemaxbutton')|strcmp(action,'zoomtimeminbutton')|strcmp(action,'zoomtimecenbutton')
    adjust_window(h.poweraxes)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoomcolbutton_update(hObject, events)
% update called from zoom buttons

action = get(hObject,'Tag');
value = str2num(get(hObject,'String'));
handles = guihandles(gcf);
F  = spm_figure('GetWin','Graphics');
h = guidata(F);

switch action
    case 'zoomcolmaxbutton'
        h.ColMax=num2str(value);
        h.ColCen=num2str(mean([str2num(h.ColMin) str2num(h.ColMax)]));
        set(handles.zoomcolcenbutton,'String',h.ColCen)
    case 'zoomcolminbutton'
        h.ColMin=num2str(value);
        h.ColCen=num2str(mean([str2num(h.ColMin) str2num(h.ColMax)]));
        set(handles.zoomcolcenbutton,'String',h.ColCen)
    case 'zoomcolcenbutton'
        h.ColMin=num2str(str2num(h.ColMin)+(value-str2num(h.ColCen)));
        set(handles.zoomcolminbutton,'String',h.ColMin)
        h.ColMax=num2str(str2num(h.ColMax)+(value-str2num(h.ColCen)));
        set(handles.zoomcolmaxbutton,'String',h.ColMax)
        h.ColCen=num2str(value);
end
guidata(F,h)

% make plots
if strcmp(action,'zoomcolmaxbutton') || strcmp(action,'zoomcolminbutton') || strcmp(action,'zoomcolcenbutton')
    draw_colors(h.poweraxes)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventbutton_update(hObject, events)
% update called from zoom buttons
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
action = get(hObject,'Value');

h.action=action;
guidata(F,h);

% make plots
draw_bars(action)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function channelslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

% Update display of current trial number
if ind>1
    SelChan=ind-1;
    set(handles.channeltext, 'String', mat2str(SelChan));
else
    SelChan=1:D.nchannels;
    set(handles.channeltext, 'String', 'All');
end

h.SelChan=SelChan;
guidata(F,h);

%Store data to plot
if get(handles.normalisebutton,'Value')
    Data=h.DataNorm;
else
    Data=h.DataDisp;
end
data_plots(Data,1)

% make plots
draw_plots(1,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trialslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

% Update display of current trial number
if ind>1
    SelTrial=ind-1;
    set(handles.trialgotext, 'String', mat2str(SelTrial));
else
    SelTrial=1:D.ntrials;
    set(handles.trialgotext, 'String', 'All');
end

h.SelTrial=SelTrial;
guidata(F,h);

%Store data to plot
if get(handles.normalisebutton,'Value')
    Data=h.DataNorm;
else
    Data=h.DataDisp;
end
data_plots(Data,1)

% Update slider for events types
Events=D.events{SelTrial};
Neventcode=max([size(Events,2) 2]);
NeventcodeStep=max([2 Neventcode]);
set(handles.eventtypeslider,'Max', Neventcode, 'Value', 1,'SliderStep', [1/(NeventcodeStep-1) min(NeventcodeStep-1, 1/(NeventcodeStep-1))]);
set(handles.eventtypetext, 'String', '1');

% make plots
draw_plots(1,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventtypeslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
if size(D,4)<1
    Events=h.Events;
else
    Events=D.events{h.SelTrial};
end
% slider value
ind = round(get(hObject, 'Value'));
tmp=size(Events,2);
if ind>tmp
    ind=tmp;
end
set(hObject, 'Value', ind);
set(hObject, 'Value', ind);

% Update display of current trial number
try
    set(handles.eventtypetext, 'String', mat2str(ind));
    set(handles.eventnametext, 'String', Events(ind).type)
end

% plot bars
try
    action = h.action;
catch
    action=1;
end
draw_bars(action)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eventslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;

ntype=str2num(get(handles.eventtypetext, 'String'));

% slider value
ind = round(get(hObject, 'Value'));
tmp=length(find(h.Code==ntype));
if ind>tmp
    ind=tmp;
end
set(hObject, 'Value', ind);

% Update display of current trial number
set(handles.eventgotext, 'String', mat2str(ind));

nevent=str2num(get(handles.eventgotext, 'String'));
tmp=find(h.Code==ntype);
Events=h.Events;
offsettime=indsample(D,Events(tmp(nevent)).time);
% set(handles.timeslider, 'Value', offset);
h.offsettime=offsettime;
h.ntype=ntype;
guidata(F,h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeslider_update(hObject, events)
% update called from zoom buttons

F = spm_figure('GetWin','Graphics');
h = guidata(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

offsettime=ind;
h.offsettime=offsettime;
guidata(F,h);

% make plots
adjust_window(h.poweraxes, h.time, h.ztime, h.zfreq, h.offsettime, h.offsetfreq)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function freqslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

offsetfreq=ind;
h.offsetfreq=offsetfreq;
guidata(F,h);

% make plots
adjust_window(h.poweraxes, h.time, h.ztime, h.zfreq, h.offsettime, h.offsetfreq)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function adjust_window(axe)
% This function adjusts data window

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

D=h.D;
time=h.time;
freq=D.tf.frequencies;
try
    TimeIndex=find(time>=str2num(h.TimeMin)&time<=str2num(h.TimeMax));
catch
    TimeIndex=1:length(time);
end
try
    FreqIndex=find(freq>=str2num(h.FreqMin)&freq<=str2num(h.FreqMax));
catch
    FreqIndex=1:length(freq);
end

axes(axe);
axis([min(time(TimeIndex)) max(time(TimeIndex)) min(freq(FreqIndex)) max(freq(FreqIndex))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_plots(FlagAxes,FlagBars)
% This function plots data in the subplots

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

D=h.D;
Events=h.Events;
SelChan=h.SelChan;

if get(handles.normalisebutton,'Value') || (min(h.DataPlot(:)) < 0)
    h.col = [flipud(fliplr(hot(128)));hot(128)];
    h.col = colormap(h.col);
else
    colormap('default');
    h.col = colormap;
end

if FlagAxes==1
    DataPlot=h.DataPlot;
    time=h.time;
    zcol=h.zcol;
    h.MaxPlot=max(abs(DataPlot(:)));
    h.MinPlot=min(abs(DataPlot(:)));
    MaxPlot=h.MaxPlot;
    MinPlot=h.MinPlot;
    axe=h.poweraxes;
elseif FlagAxes==2
    if get(handles.normalisebutton,'Value')
        DataPlot=h.DataNorm2;
    else
        DataPlot=h.DataPlot2;
    end
    time=h.time2;
    zcol=h.zcol2;
    h.MaxPlot2=max(abs(DataPlot(:)));
    h.MinPlot2=min(abs(DataPlot(:)));
    MaxPlot=h.MaxPlot2;
    MinPlot=h.MinPlot2;
    axe=h.poweraxes2;
    Av=h.Av;
end
axes(axe);
cla

imagesc(time,D.tf.frequencies,DataPlot)
axis xy
axis([min(time) max(time) min(D.tf.frequencies) max(D.tf.frequencies)])
if get(handles.normalisebutton,'Value')|| min(h.DataPlot(:))<0
    caxis([-zcol*MaxPlot zcol*MaxPlot])
else
    caxis([MinPlot zcol*MaxPlot])
end
colorbar;
xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
ylabel('Frequency [Hz]', 'Interpreter', 'tex', 'FontSize', FS1);
if FlagAxes==1
    if length(SelChan)==D.nchannels
        title([h.Flag ' averaged over channels'], 'Interpreter', 'tex', 'FontSize', FS2);
    else
        switch h.Flag
            case{'Synchrony','Coherence'}
                title([h.Flag ' between channels ' strcat(D.chanlabels{SelChan})], 'Interpreter', 'tex', 'FontSize', FS2);
            otherwise
                title([h.Flag ' on channel ' strcat(D.chanlabels{SelChan})], 'Interpreter', 'tex', 'FontSize', FS2);
        end
    end
    hold on
    if ~isempty(Events)
        M1=0;
        M2=2*max(D.tf.frequencies);
        for i3=1:h.Ntypes
            tmp=find(h.Code==i3);
            for i4=1:length(tmp)
                h.eventaxes{i3,i4}=line([Events(tmp(i4)).time Events(tmp(i4)).time],[M1 M2],[1.00001*MaxPlot 1.00001*MaxPlot],'color','m','linewidth',2,'Visible','off');
            end
        end
    end
    hold off
    guidata(F,h)
    % plot bars
    if FlagBars>0
        action = get(handles.dispeventbutton,'Value');
        if iscell(action)
            action = action{1};
        end
        draw_bars(action)
    end
elseif FlagAxes==2   
    h.eventaxes2=line([0 0],[0 2*max(D.tf.frequencies)],[1.00001*MaxPlot 1.00001*MaxPlot],'color','m','linewidth',2,'Visible','on');
    str=[h.Flag ' averaged on event(s) ' strcat(D.events.name{Av.Event}) ' and on channel(s) ' strcat(D.channels.name{Av.Chan})];
    title(str, 'Interpreter', 'tex', 'FontSize', FS2);
end

guidata(F,h)
adjust_window(axe)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_colors(hsel)
% This function plots data in the subplots

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

D=h.D;

axes(hsel);
if D.tf.rm_baseline|min(h.DataPlot(:))<0
    MaxPlot=str2num(h.ColMax);
    caxis([-MaxPlot MaxPlot])
else
    MinPlot=str2num(h.ColMin);
    MaxPlot=str2num(h.ColMax);
    caxis([MinPlot MaxPlot])
end
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_bars(action)
% This function plots events in the subplots

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
Types=h.Types;
Nevents=length(Types);
if Nevents>1
    Ntype=str2num(get(handles.eventtypetext,'String'));
else
    Ntype=1;
end

if Nevents>0
    n=length(find(h.Code==Ntype));
    for i1=1:Nevents
        set(h.eventaxes{h.Code(i1),length(find(h.Code(1:i1)==h.Code(i1)))},'Visible','off')
    end
    for i2=1:n
        if action
            set(h.eventaxes{Ntype,i2},'Visible','on')
        end
    end
end
guidata(F,h)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function averagebutton(D,SelChan)
% This function computes and display event-related average

F = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
Events=h.Events;

Fi  = spm_figure('GetWin','Interactive');
figure(Fi);clf
Av.Event=spm_input('Event type(s) ', 1, 'r',get(handles.eventtypeslider,'Value'));
tmp=get(handles.channelslider,'Value');
if tmp==1
    Av.Chan=spm_input('Channels ', '+1', 'r',1:D.nchannels);
else
    Av.Chan=spm_input('Channels ', '+1', 'r',tmp-1);
end
% Av.Freq=spm_input('Frequency window [Hz]', '+1', 'r',[D.tf.frequencies(1) D.tf.frequencies(end)]);
Av.DurPre=spm_input('Pre event duration [sec] ', '+1', 'r',10);
Av.DurPost=spm_input('Post event duration [sec] ', '+1', 'r',10);

%Compute average
Av.Onset=[];
for i1=Av.Event
    Av.Onset=[Av.Onset find(h.Code'==i1)];
end
% Av.Freq=find(D.tf.frequencies>=Av.Freq(1)&D.tf.frequencies<=Av.Freq(2));
Av.DurPre=-ceil(Av.DurPre*D.fsample);
Av.DurPost=ceil(Av.DurPost*D.fsample);
DataPlot2=zeros(D.Nfrequencies,length([Av.DurPre:Av.DurPost]));
DataOn=zeros(1,length([Av.DurPre:Av.DurPost]));
for i1=Av.Onset
    win=Events.time(i1)+[Av.DurPre:Av.DurPost];
    tmp=find(win>0&win<=D.nsamples);
    DataOn(tmp)=DataOn(tmp)+1;
    win=win(tmp);
    for i2=Av.Chan
        DataPlot2(:,tmp)=DataPlot2(:,tmp)+squeeze(h.DataDisp(i2,:,win));
    end
end
h.DataPlot2=DataPlot2./(ones(D.Nfrequencies,1)*DataOn);

%Plot average
h.time2=h.time(1-Av.DurPre+[Av.DurPre:Av.DurPost])-h.time(-Av.DurPre+1);
h.ztime2=1;
h.zfreq2=1;
h.offsettime2=find(h.time2==0);
h.offsetfreq2=1;
h.zcol2=1;
if ~isfield(h,'poweraxes2')
    h.poweraxes2 = axes('Position',...
        [0.09 0.10 0.85 0.35],...
        'Parent', F);
end
h.Av=Av;

set(handles.zoomcoloutbutton2,'Visible','on');
set(handles.zoomfreqinbutton2,'Visible','on');
set(handles.zoomtimeinbutton2,'Visible','on');

guidata(F,h)

if get(handles.normalisebutton,'Value')
    normalisebutton_update
end
draw_plots(2,1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function normalisebutton_update(varargin)
% This function computes and display frequency normalisation

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
D=h.D;
Events=h.Events;
action =get(handles.normalisebutton,'Value');
Fi  = spm_figure('GetWin','Interactive');
figure(Fi);
clf

if action
    if h.Ntypes>1
        Av.Event=spm_input('Baseline event type(s) ', 1, 'r',get(handles.eventtypeslider,'Value'));
    else
        Av.Event=1;
    end
    Av.DurPre=spm_input('Pre event duration [sec] ', '+1', 'r',10);
    Av.DurPost=spm_input('Post event duration [sec] ', '+1', 'r',0);
    Av.DurPre=-ceil(Av.DurPre*D.fsample);
    Av.DurPost=ceil(Av.DurPost*D.fsample);
    if isfield(h,'DataPlot2')
        Av.DurPre2=spm_input('Pre event duration on the mean [sec] ', '+1', 'r',10);
        Av.DurPost2=spm_input('Post event duration on the mean [sec] ', '+1', 'r',0);
        Av.DurPre2=-ceil(Av.DurPre2*D.fsample);
        Av.DurPost2=ceil(Av.DurPost2*D.fsample);
    end    

    %Compute baseline intervals
    Av.Onset=[];
    for i1=Av.Event
        Av.Onset=[Av.Onset find(h.Code'==i1)];
    end
    
    D.tf.Sbaseline = zeros(length(Av.Onset),2);
    for i1=1:length(Av.Onset)
        win=Events.time(Av.Onset(i1))+[Av.DurPre Av.DurPost];
        win(find(win<1))=1;
        win(find(win>D.nsamples))=D.nsamples;
        D.tf.Sbaseline(i1,:)=win;
    end
    D.tf.rm_baseline=1;
    Data=zeros(size(Datadisp));
    for i1=1:size(Data,1)
        Data(i1,:,:)=h.DataDisp(i1,:,:);
    end
    Data = ImaGIN_spm_eeg_bc(D, Data);
    h.DataNorm=Data;
    if isfield(h,'DataPlot2')
        tmp=find(h.time2==0);
        DD.tf.Sbaseline=tmp+[Av.DurPre2 Av.DurPost2];
        DD.tf.Sbaseline=DD.tf.Sbaseline(find(DD.tf.Sbaseline>=1&DD.tf.Sbaseline<=length(h.time2)));
        DD.Nfrequencies=D.Nfrequencies;
        DD.tf.channels=1;
        clear tmp
        tmp(1,:,:)=h.DataPlot2;
        h.DataNorm2= squeeze(ImaGIN_spm_eeg_bc(DD, tmp));
    end
else
    h=rmfield(h,'DataNorm');
    D.tf.rm_baseline=0;
    Data=h.DataDisp;
end

h.D=D;
guidata(F,h);

% make plots
data_plots(Data,1)
draw_plots(1,1)
if isfield(h,'DataPlot2')
    draw_plots(2,1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_plots(Data,Flag)
% This function computes data (average over SelChan) before plot

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);
SelChan=h.SelChan;
SelTrial=h.SelTrial;

if Flag==1
    DataPlot=zeros(size(Data,2),size(Data,3));
    for i1=SelChan
        try
            if length(SelTrial)>1
                DataPlot=DataPlot+squeeze(mean(Data(i1,:,:,SelTrial),4));
            else
                DataPlot=DataPlot+squeeze(Data(i1,:,:,SelTrial));
            end
        catch
            DataPlot=DataPlot+squeeze(Data(i1,:,:));
        end
    end
    DataPlot=DataPlot/length(SelChan);
    h.MaxPlot=max(abs(DataPlot(:)));
    h.MinPlot=min(abs(DataPlot(:)));
    h.DataPlot=DataPlot;
elseif Flag==2
    DataPlot2=zeros(size(Data,2),size(Data,3));
    for i1=SelChan
        DataPlot2=DataPlot2+squeeze(mean(Data(i1,:,:,SelTrial),4));
    end
    DataPlot2=DataPlot2/length(SelChan);
    h.MaxPlot2=max(abs(DataPlot2(:)));
    h.MinPlot2=min(abs(DataPlot2(:)));
    h.DataPlot2=DataPlot2;
end

h.ColMax=num2str(h.MaxPlot,4);
h.ColMin=num2str(h.MinPlot,4);
h.ColCen=num2str(mean([str2num(h.ColMin) str2num(h.ColMax)]),4);
try
    set(handles.zoomcolminbutton,'String',h.ColMin)
    set(handles.zoomcolcenbutton,'String',h.ColCen)
    set(handles.zoomcolmaxbutton,'String',h.ColMax)
end

guidata(F,h);

