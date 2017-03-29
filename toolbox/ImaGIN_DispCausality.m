function ImaGIN_DispCausality(varargin)
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

if nargin==0
    Fgraph  = spm_figure('GetWin','Graphics');
    Finter  = spm_figure('GetWin','Interactive');
    spm_figure('Clear',Fgraph)
    spm_figure('Clear',Finter)

    t = spm_select(inf, '\.mat$', 'Select causality data file',[],pwd,'ca1_');
    e=spm_str_manip(t(1,:),'t');
    switch e(1:2)
        case 'cc'
            ImaGIN_DispCausality_CC(t)
        case {'gs','tg'}
            ImaGIN_DispCausality_GS(t)
        case {'h2'}
            ImaGIN_DispCausality_H2(t)
        case {'ar','ta','sa'}
            ImaGIN_DispCausality_AR(t)
    end
else
    eval(varargin{1})
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_DispCausality_CC(t)

WS = spm('WinScale');

warning off

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 9);

Fgraph  = spm_figure('GetWin','Graphics');
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)


D1=spm_eeg_ldata(t);    %coupling strength

d=spm_str_manip(t,'h');
e=spm_str_manip(t,'t');
tmp=strfind(e,'_');
e=[e(1:tmp(1)+2) '2' e(tmp(1)+4:end)];
t=fullfile(d,e);
D2=spm_eeg_ldata(t);    %time lag

%plot electrodes
% tmp=spm_input('Plot electrodes ',1,'Yes|No');
% if strcmp(tmp,'Yes')
%     FlagElec=1;
Fchannels=fullfile(spm('dir'), 'EEGtemplates',D1.channels.ctf);
Template=load(Fchannels);
Pos=Template.Cpos2(:,D1.ca.channels)';
h.Pos=Pos;

%Plot glass brain
switch D1.Atlas
    case 'Human'
        P1=fullfile(spm('dir'));
    case {'Rat','Mouse','Macaque'}
        P1=fullfile(spm('dir'),'atlas5',lower(D1.Atlas));
end
path(P1,path);
% global defaults
% spm_defaults;
hMIPax = axes('Parent',Fgraph,'Position',[0.13 0.52 0.715 0.468],'Visible','off');
hMIPax = spm_mip_ui([],zeros(3,0),eye(4),[2 2 2]',hMIPax);

%-Axis offsets for 3d MIPs:
%=======================================================================
try
    Species=D1.Atlas;
catch
    Species=spm_input('Species ','+1','Human|Rat|Mouse');
end
%-MIP pane dimensions and Talairach origin offsets
%-See spm_project.c for derivation
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
switch Species
    case{'Human'}
        %             DXYZ = [182 218 182];
        %             CXYZ = [091 127 073];
        DXYZ = [192 218 182];
        CXYZ = [096 132 078];
    case{'Rat','Mouse'}
        DXYZ = [189 211 175];
        CXYZ = [106 167 141];
%         Pos=10*Pos;
end
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;

if size(D1.data,2)==1
    mas=max(D1.data(:));
    mis=min(D1.data(:));
    md=max(abs(D2.data(:)));
else
    mas=max(max(D1.data(:,:)));
    mis=min(min(D1.data(:,:)));
    md=max(max(abs(D2.data(:,:))));
end
% Zmax = str2num(spm_input(sprintf('+ Threshold synchro (max=%0.3f) ',mas), '+1', 's',0));
% Zmin = str2num(spm_input(sprintf('- Threshold synchro (min=%0.3f) ',mis), '+1', 's',0));
% Zd = str2num(spm_input(sprintf('+ Threshold delay (max=%0.3f) ',mad), '+1', 's',0));
Zmax=(mas+mis)/2;
Zmin=(mas+mis)/2;
Zd=md/2;
% if Zmax<Zmin
%     Zmin=Zmax;
% end

h.Pos=Pos;
h.Zmax=Zmax;
h.Zmin=Zmin;
h.Zd=Zd;
h.D1=D1;
h.D2=D2;
h.offset=1;
h.Po=Po;
h.hMIPax=hMIPax;


if length(D1.ca.Time)>1
    %Slider for time
    handles.timeslider=uicontrol(Fgraph, 'Tag', 'timeslider', 'Style', 'slider',...
        'Min', 1, 'Max', D1.Nsamples, 'Value', 1, 'Units',...
        'normalized', 'Position', [0.425 0.48 0.2 0.03],...
        'SliderStep', [1/(D1.Nsamples-1) unique(max([1/(D1.Nsamples-1) 0.2]))],...
        'TooltipString', 'Navigate through time',...
        'Callback', @timeslider_update,...
        'Parent', Fgraph, 'Interruptible', 'off');
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s]',...
        'Position',[0.425 0.45 0.2 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'text', 'Tag', 'timedisp',...
        'Units', 'normalized',...
        'String', num2str(D1.time(get(handles.timeslider, 'Value'))),...
        'Position',[0.425 0.511 0.2 0.029],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', spm('Colour'));
    
    %Slider for threshold (coupling, max)
    handles.couplingmaxslider=uicontrol(Fgraph, 'Tag', 'couplingmaxslider', 'Style', 'slider',...
        'Min', floor(1e3*mis), 'Max', ceil(1e3*mas), 'Value', 1e3*Zmax, 'Units',...
        'normalized', 'Position', [0.8 0.5 0.1 0.02],...
        'SliderStep', [1/(ceil(1e3*mas)-floor(1e3*mis)+1) 10/(ceil(1e3*mas)-floor(1e3*mis)+1)],...
        'TooltipString', 'Threshold applied on coupling strength (glass brain, red)',...
        'Callback', @couplingmaxslider_update,...
        'Parent', Fgraph, 'Interruptible', 'off');
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Coupling (max) ',...
        'Position',[0.65 0.5 0.15 0.02],...
        'HorizontalAlignment', 'right', 'FontSize', FS3,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'text', 'Tag', 'couplingmaxdisp',...
        'Units', 'normalized',...
        'String', num2str(1e-3*get(handles.couplingmaxslider, 'Value')),...
        'Position',[0.9 0.5 0.1 0.02],...
        'HorizontalAlignment', 'center', 'FontSize', FS3,...
        'BackgroundColor', spm('Colour'));

    %Slider for threshold (coupling, down)
    handles.couplingminslider=uicontrol(Fgraph, 'Tag', 'couplingminslider', 'Style', 'slider',...
        'Min', floor(1e3*mis), 'Max', ceil(1e3*mas), 'Value', 1e3*Zmin, 'Units',...
        'normalized', 'Position', [0.8 0.47 0.1 0.02],...
        'SliderStep', [1/(ceil(1e3*mas)-floor(1e3*mis)+1) 10/(ceil(1e3*mas)-floor(1e3*mis)+1)],...
        'TooltipString', 'Threshold applied on coupling strength (glass brain, blue)',...
        'Callback', @couplingminslider_update,...
        'Parent', Fgraph, 'Interruptible', 'off');
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Coupling (min) ',...
        'Position',[0.65 0.47 0.15 0.02],...
        'HorizontalAlignment', 'right', 'FontSize', FS3,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'text', 'Tag', 'couplingmindisp',...
        'Units', 'normalized',...
        'String', num2str(1e-3*get(handles.couplingminslider, 'Value')),...
        'Position',[0.9 0.47 0.1 0.02],...
        'HorizontalAlignment', 'center', 'FontSize', FS3,...
        'BackgroundColor', spm('Colour'));

    %Slider for threshold (delay)
    handles.delayslider=uicontrol(Fgraph, 'Tag', 'delayslider', 'Style', 'slider',...
        'Min', 0, 'Max', ceil(1e3*md), 'Value', 1e3*Zd, 'Units',...
        'normalized', 'Position', [0.8 0.44 0.1 0.02],...
        'SliderStep', [1/(ceil(1e3*md)-1) 2/(ceil(1e3*md)-1)],...
        'TooltipString', 'Threshold applied on delay (glass brain, arrow)',...
        'Callback', @delayslider_update,...
        'Parent', Fgraph, 'Interruptible', 'off');
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Delay ',...
        'Position',[0.65 0.44 0.15 0.02],...
        'HorizontalAlignment', 'right', 'FontSize', FS3,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'text', 'Tag', 'delaydisp',...
        'Units', 'normalized',...
        'String', num2str(1e-3*get(handles.delayslider, 'Value')),...
        'Position',[0.9 0.44 0.1 0.02],...
        'HorizontalAlignment', 'center', 'FontSize', FS3,...
        'BackgroundColor', spm('Colour'));

    %Axes for time series    
    for i1=1:D1.Nchannels+1
        h.Heegaxes1(i1) = axes('Position',...
            [0.13 0.26 0.83 0.15],...
            'Parent', Fgraph);
        axes(h.Heegaxes1(i1));
        if i1<=D1.Nchannels
            h.Heegplot1{i1}=plot(D1.time,D1.data(i1,:),'k');
            ylabel([sprintf('%d: ',i1) D1.ca.chanlabels{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
        else
            h.Heegplot1{i1}=plot(D1.time,D1.data(:,:));
            hold on
            h.Heegplot1{i1}=[h.Heegplot1{i1};plot(D1.time,mean(D1.data(:,:)),'LineWidth',4,'Color','k')];
            hold off
            ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
        end
        axis tight
        title('Coupling strength','FontSize', FS1)
%         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
        set(h.Heegaxes1(i1),'Visible','off')
        set(h.Heegplot1{i1},'Visible','off')
        h.eventaxes1{i1,1}=line([D1.time(1) D1.time(end)],[Zmax Zmax],'color','k','linewidth',1,'linestyle',':','Visible','off');
        h.eventaxes1{i1,2}=line([D1.time(1) D1.time(end)],[Zmin Zmin],'color','k','linewidth',1,'linestyle',':','Visible','off');
        h.eventaxes1{i1,3}=line([D1.time(h.offset) D1.time(h.offset)],[mis mas],'color','k','linewidth',1,'Visible','off');
    end
    for i1=1:D2.Nchannels+1
        h.Heegaxes2(i1) = axes('Position',...
            [0.13 0.05 0.83 0.15],...
            'Parent', Fgraph);
        axes(h.Heegaxes2(i1));
        if i1<=D2.Nchannels
            h.Heegplot2{i1}=plot(D2.time,D2.data(i1,:),'k');
            ylabel([sprintf('%d: ',i1) D2.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
        else
            h.Heegplot2{i1}=plot(D2.time,D2.data(:,:));
            hold on
            h.Heegplot2{i1}=[h.Heegplot2{i1};plot(D2.time,mean(D2.data(:,:)),'LineWidth',4,'Color','k')];
            hold off
            ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
        end
        axis tight
        title('Delay','FontSize', FS1)
        xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
        set(h.Heegaxes2(i1),'Visible','off')
        set(h.Heegplot2{i1},'Visible','off')
        h.eventaxes2{i1,1}=line([D2.time(1) D2.time(end)],[Zd Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
        h.eventaxes2{i1,2}=line([D2.time(1) D2.time(end)],[-Zd -Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
        h.eventaxes2{i1,3}=line([D2.time(h.offset) D2.time(h.offset)],[-md md],'color','k','linewidth',1,'Visible','off');
    end

    txt_box = D1.ca.chanlabels;
    txt_box{end+1} = 'all';
    h.plotelec = uicontrol(Fgraph,'Style','popup','String',txt_box, ...
        'Position',[50 400 100 20].*WS, 'Value',length(txt_box),...
        'Callback',@draw_axes);
        
    guidata(Fgraph,h);
    draw_axes;    
        
end



colordot='k';
axes(hMIPax)
for i1=1:size(Pos,1)
    hold on
    h1=plot(Po(1) + Pos(i1,2),Po(2) + Pos(i1,1),'.');
    h2=plot(Po(1) + Pos(i1,2),Po(3) - Pos(i1,3),'.');
    h3=plot(Po(4) + Pos(i1,1),Po(3) - Pos(i1,3),'.');
    hold off
    set(h1,'MarkerSize',15,'color',colordot);
    set(h2,'MarkerSize',15,'color',colordot);
    set(h3,'MarkerSize',15,'color',colordot);
end
% hquiv=zeros(6,size(D1.ca.CM,2));
% hold on
% for i1=1:size(hquiv,2)
%     hquiv(1,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(2,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(3,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(4,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(5,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(6,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
% end
% hold off

% h.hquiv=hquiv;
% store handles
guidata(Fgraph, h);

draw_lines;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_DispCausality_AR(t)

WS = spm('WinScale');

warning off

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 9);

Fgraph  = spm_figure('GetWin','Graphics');
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)

if size(t,1)==1
    clear D1
    D1=spm_eeg_load(t);
    try
        D1.ca.S(find(isnan(D1.ca.S)))=0;
    end
    try
        D1.ca.ffDTF(find(isnan(D1.ca.ffDTF)))=0;
    end
    try
        D1.ca.GdDTF(find(isnan(D1.ca.GdDTF)))=0;
    end
    try
        D1.ca.ffPDC(find(isnan(D1.ca.ffPDC)))=0;
    end
    try
        D1.ca.GffPDC(find(isnan(D1.ca.GffPDC)))=0;
    end
    D1.ca.DTF(find(isnan(D1.ca.DTF)))=0;
    D1.ca.dDTF(find(isnan(D1.ca.dDTF)))=0;
    D1.ca.PDC(find(isnan(D1.ca.PDC)))=0;
    D1.ca.GPDC(find(isnan(D1.ca.GPDC)))=0;
else
    clear D1 D2
    for i1=1:size(t,1)
        D2{i1}=spm_eeg_ldata(deblank(t(i1,:)));
%         D2{i1}.ca=rmfield(D2{i1}.ca,'ffDTF');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'GdDTF');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'ffPDC');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'GffPDC');
        D2{i1}.ca=rmfield(D2{i1}.ca,'DTF');
        D2{i1}.ca=rmfield(D2{i1}.ca,'PDC');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'ffDTF_S');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'GdDTF_S');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'ffPDC_S');
%         D2{i1}.ca=rmfield(D2{i1}.ca,'GffPDC_S');
        D2{i1}.ca=rmfield(D2{i1}.ca,'DTF_S');
        D2{i1}.ca=rmfield(D2{i1}.ca,'PDC_S');
%         D2{i1}.ca.S(find(isnan(D2{i1}.ca.S)))=0;
%         D2{i1}.ca.ffDTF(find(isnan(D2{i1}.ca.ffDTF)))=0;
%         D2{i1}.ca.GdDTF(find(isnan(D2{i1}.ca.GdDTF)))=0;
%         D2{i1}.ca.ffPDC(find(isnan(D2{i1}.ca.ffPDC)))=0;
%         D2{i1}.ca.GffPDC(find(isnan(D2{i1}.ca.GffPDC)))=0;
%         D2{i1}.ca.DTF(find(isnan(D2{i1}.ca.DTF)))=0;
        D2{i1}.ca.dDTF(find(isnan(D2{i1}.ca.dDTF)))=0;
%         D2{i1}.ca.PDC(find(isnan(D2{i1}.ca.PDC)))=0;
        D2{i1}.ca.GPDC(find(isnan(D2{i1}.ca.GPDC)))=0;
    end
    D1=D2{1};
    D1=rmfield(D1,'data');
    D1.data=D2{1}.data(:,:);
%     figure
%     subplot(211)
%     plot(D2{1}.data(1,:)')
%     subplot(212)
%     plot(D1.data(1,:)')
    for i1=2:size(t,1)
        try
            D1.data=D1.data+D2{i1}.data(:,:);
        end
        try
            D1.ca.S=D1.ca.S+D2{i1}.ca.S;
%             D1.ca.ffDTF=D1.ca.ffDTF+D2{i1}.ca.ffDTF;
%             D1.ca.GdDTF=D1.ca.GdDTF+D2{i1}.ca.GdDTF;
%             D1.ca.ffPDC=D1.ca.ffPDC+D2{i1}.ca.ffPDC;
%             D1.ca.GffPDC=D1.ca.GffPDC+D2{i1}.ca.GffPDC;
        end
%         D1.ca.DTF=D1.ca.DTF+D2{i1}.ca.DTF;
        D1.ca.dDTF=D1.ca.dDTF+D2{i1}.ca.dDTF;
%         D1.ca.PDC=D1.ca.PDC+D2{i1}.ca.PDC;
        D1.ca.GPDC=D1.ca.GPDC+D2{i1}.ca.GPDC;
        try
            D1.ca.S_S=D1.ca.S_S+D2{i1}.ca.S_S;
%             D1.ca.ffDTF_S=D1.ca.ffDTF_S+D2{i1}.ca.ffDTF_S;
%             D1.ca.GdDTF_S=D1.ca.GdDTF_S+D2{i1}.ca.GdDTF_S;
%             D1.ca.ffPDC_S=D1.ca.ffPDC_S+D2{i1}.ca.ffPDC_S;
%             D1.ca.GffPDC_S=D1.ca.GffPDC_S+D2{i1}.ca.GffPDC_S;
        end
%         D1.ca.DTF_S=D1.ca.DTF_S+D2{i1}.ca.DTF_S;
        D1.ca.dDTF_S=D1.ca.dDTF_S+D2{i1}.ca.dDTF_S;
%         D1.ca.PDC_S=D1.ca.PDC_S+D2{i1}.ca.PDC_S;
        D1.ca.GPDC_S=D1.ca.GPDC_S+D2{i1}.ca.GPDC_S;
%         figure
%         subplot(211)
%         plot(D2{i1}.data(1,:)')
%         subplot(212)
%         plot(D1.data(1,:)')
    end
    D1.data(:,:)=D1.data(:,:)./size(t,1);
    try
        D1.ca.S=D1.ca.S./size(t,1);
%         D1.ca.ffDTF=D1.ca.ffDTF./size(t,1);
%         D1.ca.GdDTF=D1.ca.GdDTF./size(t,1);
%         D1.ca.ffPDC=D1.ca.ffPDC./size(t,1);
%         D1.ca.GffPDC=D1.ca.GffPDC./size(t,1);
    end
%     D1.ca.DTF=D1.ca.DTF./size(t,1);
    D1.ca.dDTF=D1.ca.dDTF./size(t,1);
%     D1.ca.PDC=D1.ca.PDC./size(t,1);
    D1.ca.GPDC=D1.ca.GPDC./size(t,1);
    try
        D1.ca.S_S=D1.ca.S_S./size(t,1);
%         D1.ca.ffDTF_S=D1.ca.ffDTF_S./size(t,1);
%         D1.ca.GdDTF_S=D1.ca.GdDTF_S./size(t,1);
%         D1.ca.ffPDC_S=D1.ca.ffPDC_S./size(t,1);
%         D1.ca.GffPDC_S=D1.ca.GffPDC_S./size(t,1);
    end
%     D1.ca.DTF_S=D1.ca.DTF_S./size(t,1);
    D1.ca.dDTF_S=D1.ca.dDTF_S./size(t,1);
%     D1.ca.PDC_S=D1.ca.PDC_S./size(t,1);
    D1.ca.GPDC_S=D1.ca.GPDC_S./size(t,1);
end

DataPlot1=squeeze(mean(D1.ca.dDTF,4));
DataPlot1=permute(DataPlot1,[2 1 3]);
DataPlot2=DataPlot1-permute(DataPlot1,[2 1 3]);       
h.DataPlot1=DataPlot1;
h.DataPlot2=DataPlot2;
h.Mode='Time';
h.RangeTime=[min(D1.ca.Time) max(D1.ca.Time)];
h.RangeTimeIndex=1:length(D1.ca.Time);
h.RangeFreq=[min(D1.ca.Freq) max(D1.ca.Freq)];
h.RangeFreqIndex=1:length(D1.ca.Freq);
h.DataType='dDTF';

%plot electrodes
tmp=sensors(D1,'EEG');
Pos=tmp.pnt(D1.ca.channels,:);
h.Pos=Pos;

%Plot glass brain
if ~isfield(D1,'Atlas')
    D1.Atlas='Human';
end
switch D1.Atlas
    case 'Human'
        P1=spm('dir');
    case {'Rat','Mouse','Macaque'}
        P1=fullfile(spm('dir'),'atlas5',lower(D1.Atlas));
end
path(P1,path);
% global defaults
% spm_defaults;
spm_figure('Clear',Fgraph)
hMIPax = axes('Parent',Fgraph,'Position',[0.13 0.52 0.715 0.468],'Visible','off');
hMIPax = spm_mip_ui([],zeros(3,0),eye(4),[2 2 2]',hMIPax);
h.hMIPax=hMIPax;

%-Axis offsets for 3d MIPs:
%=======================================================================
try
    Species=D1.Atlas;
catch
    Species=spm_input('Species ','+1','Human|Rat|Mouse');
end
%-MIP pane dimensions and Talairach origin offsets
%-See spm_project.c for derivation
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
switch Species
    case{'Human'}
        %             DXYZ = [182 218 182];
        %             CXYZ = [091 127 073];
        DXYZ = [192 218 182];
        CXYZ = [096 132 078];
    case{'Rat','Mouse'}
        DXYZ = [189 211 175];
        CXYZ = [106 167 141];
%         Pos=10*Pos;
end
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;
h.Po=Po;

if nsamples(D1)==1
    mas=max(D1(:,:));
    mis=min(D1(:,:));
else
    mas=max(max(D1(:,:)));
    mis=min(min(D1(:,:)));
end
% Zmax=(mas+mis)/2;
% Zmin=(mas+mis)/2;
% Zd=md/2;
% h.Zmax=Zmax;
% h.Zmin=Zmin;
% h.Zd=Zd;

h.D1=D1;
h.offset1=1;
h.offset2=1;
h.pvalue=min(D1.ca.p)+eps;

if length(D1.ca.Time)>1
%     %Slider for time
%     handles.timeslider=uicontrol(Fgraph, 'Tag', 'timeslider', 'Style', 'slider',...
%         'Min', 1, 'Max', D1.Nsamples, 'Value', 1, 'Units',...
%         'normalized', 'Position', [0.425 0.48 0.2 0.03],...
%         'SliderStep', [1/(D1.Nsamples-1) unique(max([1/(D1.Nsamples-1) 0.2]))],...
%         'TooltipString', 'Navigate through time',...
%         'Callback', @timeslider_update,...
%         'Parent', Fgraph, 'Interruptible', 'off');
%     uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
%         'String', 'Time [s]',...
%         'Position',[0.425 0.45 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', 'w');
%     uicontrol(Fgraph, 'Style', 'text', 'Tag', 'timedisp',...
%         'Units', 'normalized',...
%         'String', num2str(D1.time(get(handles.timeslider, 'Value'))),...
%         'Position',[0.425 0.511 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', spm('Colour'));
    
    %p value update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'p value ',...
        'Position',[0.2 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'pvalue',...
        'Units', 'normalized',...
        'String', num2str(h.pvalue),...
        'Position',[0.3 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @pvalue_update,...
        'BackgroundColor', spm('Colour'));
    
    %range update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] ',...
        'Position',[0.58 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'rangetime',...
        'Units', 'normalized',...
        'String', num2str(h.RangeTime),...
        'Position',[0.68 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @rangetime_update,...
        'BackgroundColor', spm('Colour'));
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Freq [Hz] ',...
        'Position',[0.79 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'rangefreq',...
        'Units', 'normalized',...
        'String', num2str(h.RangeFreq),...
        'Position',[0.89 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @rangefreq_update,...
        'BackgroundColor', spm('Colour'));
    
    %Display Time
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] = ',...
        'Position',[0.55 0.55 0.15 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'currenttime',...
        'Units', 'normalized',...
        'String', num2str(h.D1.ca.Time(h.offset1)),...
        'Position',[0.69 0.55 0.08 0.028],...
        'HorizontalAlignment', 'left', 'FontSize', FS2,...
        'Callback', @currenttime_update,...
        'BackgroundColor', spm('Colour'));

    %Movie
    uicontrol(Fgraph, 'Style', 'pushbutton', 'Tag', 'movie',...
        'Units', 'normalized',...
        'String', 'Movie',...
        'Position',[0.9 0.55 0.08 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @makemovie,...
        'BackgroundColor', spm('Colour'));
    

    %Axes for ar data
    SizeAxeX=0.55/nchannels(D1);
    SizeAxeY=0.43/nchannels(D1);
    for i1=1:nchannels(D1)
        for i2=1:nchannels(D1)
            h.Heegaxes1(i1,i2) = axes('Position',[0.07+(i1-1)*SizeAxeX 0.02+(nchannels(D1)-i2)*SizeAxeY 0.95*SizeAxeX 0.95*SizeAxeY],...
                'Parent', Fgraph);
%             axes(h.Heegaxes1(i1,i2));
%             area(D1.ca.Time,squeeze(DataPlot1(i1,i2,:)));
%             h.eventaxes1(i1,i2)=line([D1.time(h.offset1) D1.time(h.offset1)],[0 1],'color','m','linewidth',1,'Visible','off');
        end
    end
    h.SizeAxeX=SizeAxeX;
    h.SizeAxeY=SizeAxeY;

    
    %Axes for zoom
    h.Heegaxes2(1) = axes('Position',[0.66 0.34 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(1));
%     h.eventaxes2{1} = line([D1.time(h.offset2) D1.time(h.offset2)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(2) = axes('Position',[0.66 0.19 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(2));
%     h.eventaxes2{2} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(3) = axes('Position',[0.66 0.04 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(3));
%     h.eventaxes2{3} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    
% 
%         if i1<=nchannels(D1)
%             h.Heegplot1{i1}=plot(D1.time,D1.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D1.ca.chanlabels{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot1{i1}=plot(D1.time,D1.data(:,:));
%             hold on
%             h.Heegplot1{i1}=[h.Heegplot1{i1};plot(D1.time,mean(D1.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Coupling strength','FontSize', FS1)
% %         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes1(i1),'Visible','off')
%         set(h.Heegplot1{i1},'Visible','off')
%         h.eventaxes1{i1,1}=line([D1.time(1) D1.time(end)],[Zmax Zmax],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,2}=line([D1.time(1) D1.time(end)],[Zmin Zmin],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,3}=line([D1.time(h.offset) D1.time(h.offset)],[mis mas],'color','k','linewidth',1,'Visible','off');
%     end

%     for i1=1:D2.Nchannels+1
%         h.Heegaxes2(i1) = axes('Position',...
%             [0.13 0.05 0.83 0.15],...
%             'Parent', Fgraph);
%         axes(h.Heegaxes2(i1));
%         if i1<=D2.Nchannels
%             h.Heegplot2{i1}=plot(D2.time,D2.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D2.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot2{i1}=plot(D2.time,D2.data(:,:));
%             hold on
%             h.Heegplot2{i1}=[h.Heegplot2{i1};plot(D2.time,mean(D2.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Delay','FontSize', FS1)
%         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes2(i1),'Visible','off')
%         set(h.Heegplot2{i1},'Visible','off')
%         h.eventaxes2{i1,1}=line([D2.time(1) D2.time(end)],[Zd Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,2}=line([D2.time(1) D2.time(end)],[-Zd -Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,3}=line([D2.time(h.offset) D2.time(h.offset)],[-md md],'color','k','linewidth',1,'Visible','off');
%     end

    txt_box{1}='Time';
    txt_box{2} = 'Freq';
    h.change_mode = uicontrol(Fgraph,'Style','popup','String',txt_box, ...
        'Position',[280 430 70 20].*WS, 'Value',1,...
        'Callback',@change_mode);
        
    txt_box_data{1} = 'PSD';
    txt_box_data{2} = 'DTF';
    txt_box_data{3} = 'ffDTF';
    txt_box_data{4} = 'dDTF';
    txt_box_data{5} = 'GdTF';
    txt_box_data{6} = 'PDC';
    txt_box_data{7} = 'ffPDC';
    txt_box_data{8} = 'GPDC';
    txt_box_data{9} = 'GffPDC';
    h.load_data = uicontrol(Fgraph,'Style','popup','String',txt_box_data, ...
        'Position',[50 430 70 20].*WS, 'Value',4,...
        'Callback',@load_data);

        
    h.PairSelect=[1 1];
    guidata(Fgraph,h);
    [h.Threshold1,h.Threshold2]=pvalue_compute;
    guidata(Fgraph,h);
    h=draw_axes_ar;    
        
end


colordot='k';
axes(hMIPax)
for i1=1:size(Pos,1)
    hold on
    h1=plot(Po(1) + Pos(i1,2),Po(2) + Pos(i1,1),'.');
    h2=plot(Po(1) + Pos(i1,2),Po(3) - Pos(i1,3),'.');
    h3=plot(Po(4) + Pos(i1,1),Po(3) - Pos(i1,3),'.');
    hold off
    set(h1,'MarkerSize',15,'color',colordot);
    set(h2,'MarkerSize',15,'color',colordot);
    set(h3,'MarkerSize',15,'color',colordot);
end
% hquiv=zeros(6,size(D1.ca.CM,2));
% hold on
% for i1=1:size(hquiv,2)
%     hquiv(1,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(2,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(3,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(4,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(5,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(6,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
% end
% hold off

% h.hquiv=hquiv;
% store handles
guidata(Fgraph, h);

% draw_lines;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_DispCausality_H2(t)

WS = spm('WinScale');

warning off

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 9);

Fgraph  = spm_figure('GetWin','Graphics');
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)

if size(t,1)==1
    clear D1
    D1=spm_eeg_load(t);
    D1.ca.H2X(find(isnan(D1.ca.H2X)))=0;
    D1.ca.H2Y(find(isnan(D1.ca.H2Y)))=0;
    D1.ca.H2=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.Time));
    D1.ca.p=[1:size(D1.ca.H2X_S,2)]/(1+size(D1.ca.H2X_S,2));
    D1.ca.H2_S=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.p));
    for i2=1:size(D1.ca.CM,2)
        D1.ca.H2(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.H2X(i2,:);
        D1.ca.H2(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.H2Y(i2,:);
        D1.ca.H2_S(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.H2X_S(i2,:);
        D1.ca.H2_S(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.H2Y_S(i2,:);
    end
else
    clear D1 D2
    for i1=1:size(t,1)
        D2{i1}=spm_eeg_load(deblank(t(i1,:)));
        D2{i1}.ca.H2X(find(isnan(D2{i1}.ca.H2X)))=0;
        D2{i1}.ca.H2Y(find(isnan(D2{i1}.ca.H2Y)))=0;
    end
    D1=D2{1};
    D1.ca.p=[1:size(D1.ca.H2X_S,2)]/(1+size(D1.ca.H2X_S,2));
    D1.ca.H2=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.Time));
    D1.ca.H2_S=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.p));
    for i1=1:size(t,1)
        for i2=1:size(D1.ca.CM,2)
            D1.ca.H2(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.H2X(i2,:)+D2{i1}.ca.H2X(i2,:);
            D1.ca.H2(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.H2Y(i2,:)+D2{i1}.ca.H2Y(i2,:);
            D1.ca.H2_S(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.H2X_S(i2,:)+D2{i1}.ca.H2X_S(i2,:);
            D1.ca.H2_S(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.H2Y_S(i2,:)+D2{i1}.ca.H2Y_S(i2,:);
        end
    end
    D1.ca.H2=D1.ca.H2./size(t,1);
    D1.ca.H2_S=D1.ca.H2_S./size(t,1);
end

DataPlot1=D1.ca.H2;
DataPlot1=permute(DataPlot1,[2 1 3]);
DataPlot2=DataPlot1-permute(DataPlot1,[2 1 3]);       
h.DataPlot1=DataPlot1;
h.DataPlot2=DataPlot2;
h.RangeTime=[min(D1.ca.Time) max(D1.ca.Time)];
h.RangeTimeIndex=1:length(D1.ca.Time);
h.Mode='Time';

%plot electrodes
% Fchannels=fullfile(spm('dir'), 'EEGtemplates',D1.channels.ctf);
% Template=load(Fchannels);
% Pos=Template.Cpos2(:,D1.ca.channels)';
TMP=sensors(D1,'EEG');
Pos=TMP.pnt(D1.ca.channels,:);
h.Pos=Pos;

%rename electrodes
% Name=cell(size(D1.channels.name));
% Name=tmp.label;
% ok=0;
% for i1=1:length(D1.channels.name)
%     tmp=findstr(D1.channels.name{i1},'/');
%     if ~isempty(tmp)
%         ok=1;
%         Name{D1.ca.CM(1,i1)}=D1.channels.name{i1}(1:tmp-1);
%         Name{D1.ca.CM(2,i1)}=D1.channels.name{i1}(tmp+1:end);
%     end
% end
% if ok
%     D1.channels.name=Name;
% end
        

%Plot glass brain
D1.Atlas='Rat';
switch D1.Atlas
    case 'Human'
        P1=spm('dir');
    case {'Rat','Mouse','Macaque'}
        P1=fullfile(spm('dir'),'atlas8',lower(D1.Atlas));
end
path(P1,path);
% global defaults
% spm_defaults;
spm_figure('Clear',Fgraph)
hMIPax = axes('Parent',Fgraph,'Position',[0.13 0.52 0.715 0.468],'Visible','off');
hMIPax = spm_mip_ui([],zeros(3,0),eye(4),[2 2 2]',hMIPax);
h.hMIPax=hMIPax;

%-Axis offsets for 3d MIPs:
%=======================================================================
try
    Species=D1.Atlas;
catch
    Species=spm_input('Species ','+1','Human|Rat|Mouse');
end
%-MIP pane dimensions and Talairach origin offsets
%-See spm_project.c for derivation
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
switch Species
    case{'Human'}
        %             DXYZ = [182 218 182];
        %             CXYZ = [091 127 073];
        DXYZ = [192 218 182];
        CXYZ = [096 132 078];
    case{'Rat','Mouse'}
        DXYZ = [189 211 175];
        CXYZ = [106 167 141];
%         Pos=10*Pos;
end
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;
h.Po=Po;

h.D1=D1;
h.offset1=1;
h.offset2=1;
h.pvalue=min(D1.ca.p)+eps;

if length(D1.ca.Time)>1
%     %Slider for time
%     handles.timeslider=uicontrol(Fgraph, 'Tag', 'timeslider', 'Style', 'slider',...
%         'Min', 1, 'Max', D1.Nsamples, 'Value', 1, 'Units',...
%         'normalized', 'Position', [0.425 0.48 0.2 0.03],...
%         'SliderStep', [1/(D1.Nsamples-1) unique(max([1/(D1.Nsamples-1) 0.2]))],...
%         'TooltipString', 'Navigate through time',...
%         'Callback', @timeslider_update,...
%         'Parent', Fgraph, 'Interruptible', 'off');
%     uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
%         'String', 'Time [s]',...
%         'Position',[0.425 0.45 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', 'w');
%     uicontrol(Fgraph, 'Style', 'text', 'Tag', 'timedisp',...
%         'Units', 'normalized',...
%         'String', num2str(D1.time(get(handles.timeslider, 'Value'))),...
%         'Position',[0.425 0.511 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', spm('Colour'));
    
    %p value update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'p value ',...
        'Position',[0.2 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'pvalue',...
        'Units', 'normalized',...
        'String', num2str(h.pvalue),...
        'Position',[0.3 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @pvalue_update,...
        'BackgroundColor', spm('Colour'));
    
    %range update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] ',...
        'Position',[0.58 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'rangetime',...
        'Units', 'normalized',...
        'String', num2str(h.RangeTime),...
        'Position',[0.68 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @rangetime_update,...
        'BackgroundColor', spm('Colour'));
    
    %Display Time
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] = ',...
        'Position',[0.55 0.55 0.15 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'currenttime',...
        'Units', 'normalized',...
        'String', num2str(h.D1.ca.Time(h.offset1)),...
        'Position',[0.69 0.55 0.08 0.028],...
        'HorizontalAlignment', 'left', 'FontSize', FS2,...
        'Callback', @currenttime_update,...
        'BackgroundColor', spm('Colour'));

    %Movie
    uicontrol(Fgraph, 'Style', 'pushbutton', 'Tag', 'movie',...
        'Units', 'normalized',...
        'String', 'Movie',...
        'Position',[0.9 0.55 0.08 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @makemovie,...
        'BackgroundColor', spm('Colour'));
    

    %Axes for ar data
    SizeAxeX=0.55/length(D1.ca.channels);
    SizeAxeY=0.43/length(D1.ca.channels);
    for i1=1:length(D1.ca.channels)
        for i2=1:length(D1.ca.channels)
            h.Heegaxes1(i1,i2) = axes('Position',[0.07+(i1-1)*SizeAxeX 0.02+(length(D1.ca.channels)-i2)*SizeAxeY 0.95*SizeAxeX 0.95*SizeAxeY],...
                'Parent', Fgraph);
%             axes(h.Heegaxes1(i1,i2));
%             area(D1.ca.Time,squeeze(DataPlot1(i1,i2,:)));
%             h.eventaxes1(i1,i2)=line([D1.time(h.offset1) D1.time(h.offset1)],[0 1],'color','m','linewidth',1,'Visible','off');
        end
    end
    h.SizeAxeX=SizeAxeX;
    h.SizeAxeY=SizeAxeY;

    
    %Axes for zoom
    h.Heegaxes2(1) = axes('Position',[0.66 0.34 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(1));
%     h.eventaxes2{1} = line([D1.time(h.offset2) D1.time(h.offset2)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(2) = axes('Position',[0.66 0.19 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(2));
%     h.eventaxes2{2} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(3) = axes('Position',[0.66 0.04 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(3));
%     h.eventaxes2{3} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    
% 
%         if i1<=nchannels(D1)
%             h.Heegplot1{i1}=plot(D1.time,D1.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D1.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot1{i1}=plot(D1.time,D1.data(:,:));
%             hold on
%             h.Heegplot1{i1}=[h.Heegplot1{i1};plot(D1.time,mean(D1.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Coupling strength','FontSize', FS1)
% %         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes1(i1),'Visible','off')
%         set(h.Heegplot1{i1},'Visible','off')
%         h.eventaxes1{i1,1}=line([D1.time(1) D1.time(end)],[Zmax Zmax],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,2}=line([D1.time(1) D1.time(end)],[Zmin Zmin],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,3}=line([D1.time(h.offset) D1.time(h.offset)],[mis mas],'color','k','linewidth',1,'Visible','off');
%     end

%     for i1=1:D2.Nchannels+1
%         h.Heegaxes2(i1) = axes('Position',...
%             [0.13 0.05 0.83 0.15],...
%             'Parent', Fgraph);
%         axes(h.Heegaxes2(i1));
%         if i1<=D2.Nchannels
%             h.Heegplot2{i1}=plot(D2.time,D2.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D2.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot2{i1}=plot(D2.time,D2.data(:,:));
%             hold on
%             h.Heegplot2{i1}=[h.Heegplot2{i1};plot(D2.time,mean(D2.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Delay','FontSize', FS1)
%         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes2(i1),'Visible','off')
%         set(h.Heegplot2{i1},'Visible','off')
%         h.eventaxes2{i1,1}=line([D2.time(1) D2.time(end)],[Zd Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,2}=line([D2.time(1) D2.time(end)],[-Zd -Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,3}=line([D2.time(h.offset) D2.time(h.offset)],[-md md],'color','k','linewidth',1,'Visible','off');
%     end
                
    h.PairSelect=[1 1];
    guidata(Fgraph,h);
    [h.Threshold1,h.Threshold2]=pvalue_compute;
    guidata(Fgraph,h);
    h=draw_axes_gs;    
        
end


colordot='k';
axes(hMIPax)
for i1=1:size(Pos,1)
    hold on
    h1=plot(Po(1) + Pos(i1,2),Po(2) + Pos(i1,1),'.');
    h2=plot(Po(1) + Pos(i1,2),Po(3) - Pos(i1,3),'.');
    h3=plot(Po(4) + Pos(i1,1),Po(3) - Pos(i1,3),'.');
    hold off
    set(h1,'MarkerSize',15,'color',colordot);
    set(h2,'MarkerSize',15,'color',colordot);
    set(h3,'MarkerSize',15,'color',colordot);
end
% hquiv=zeros(6,size(D1.ca.CM,2));
% hold on
% for i1=1:size(hquiv,2)
%     hquiv(1,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(2,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(3,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(4,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(5,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(6,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
% end
% hold off

% h.hquiv=hquiv;
% store handles
guidata(Fgraph, h);

% draw_lines;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImaGIN_DispCausality_GS(t)

WS = spm('WinScale');

warning off

FS1 = spm('FontSize', 14);
FS2 = spm('FontSize', 12);
FS3 = spm('FontSize', 9);

Fgraph  = spm_figure('GetWin','Graphics');
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Fgraph)
spm_figure('Clear',Finter)

if size(t,1)==1
    clear D1
    D1=spm_eeg_ldata(t);
    D1.ca.GSX(find(isnan(D1.ca.GSX)))=0;
    D1.ca.GSY(find(isnan(D1.ca.GSY)))=0;
    D1.ca.GS=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.Time));
    D1.ca.p=[1:size(D1.ca.GSX_S,2)]/(1+size(D1.ca.GSX_S,2));
    D1.ca.GS_S=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.p));
    for i2=1:size(D1.ca.CM,2)
        D1.ca.GS(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.GSX(i2,:);
        D1.ca.GS(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.GSY(i2,:);
        D1.ca.GS_S(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.GSX_S(i2,:);
        D1.ca.GS_S(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.GSY_S(i2,:);
    end
else
    clear D1 D2
    for i1=1:size(t,1)
        D2{i1}=spm_eeg_ldata(deblank(t(i1,:)));
        D2{i1}.ca.GSX(find(isnan(D2{i1}.ca.GSX)))=0;
        D2{i1}.ca.GSY(find(isnan(D2{i1}.ca.GSY)))=0;
    end
    D1=D2{1};
    D1.ca.p=[1:size(D1.ca.GSX_S,2)]/(1+size(D1.ca.GSX_S,2));
    D1.ca.GS=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.Time));
    D1.ca.GS_S=zeros(length(D1.ca.channels),length(D1.ca.channels),length(D1.ca.p));
    for i1=1:size(t,1)
        for i2=1:size(D1.ca.CM,2)
            D1.ca.GS(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.GSX(i2,:)+D2{i1}.ca.GSX(i2,:);
            D1.ca.GS(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.GSY(i2,:)+D2{i1}.ca.GSY(i2,:);
            D1.ca.GS_S(D1.ca.CM(1,i2),D1.ca.CM(2,i2),:)=D1.ca.GSX_S(i2,:)+D2{i1}.ca.GSX_S(i2,:);
            D1.ca.GS_S(D1.ca.CM(2,i2),D1.ca.CM(1,i2),:)=D1.ca.GSY_S(i2,:)+D2{i1}.ca.GSY_S(i2,:);
        end
    end
    D1.ca.GS=D1.ca.GS./size(t,1);
    D1.ca.GS_S=D1.ca.GS_S./size(t,1);
end

DataPlot1=D1.ca.GS;
DataPlot1=permute(DataPlot1,[2 1 3]);
DataPlot2=DataPlot1-permute(DataPlot1,[2 1 3]);       
h.DataPlot1=DataPlot1;
h.DataPlot2=DataPlot2;
h.RangeTime=[min(D1.ca.Time) max(D1.ca.Time)];
h.RangeTimeIndex=1:length(D1.ca.Time);
h.Mode='Time';

%plot electrodes
Fchannels=fullfile(spm('dir'), 'EEGtemplates',D1.channels.ctf);
Template=load(Fchannels);
Pos=Template.Cpos2(:,D1.ca.channels)';
h.Pos=Pos;

% %rename electrodes
% Name=cell(size(D1.channels.name));
% ok=0;
% for i1=1:length(D1.channels.name)
%     tmp=findstr(D1.channels.name{i1},'/');
%     if ~isempty(tmp)
%         ok=1;
%         Name{D1.ca.CM(1,i1)}=D1.channels.name{i1}(1:tmp-1);
%         Name{D1.ca.CM(2,i1)}=D1.channels.name{i1}(tmp+1:end);
%     end
% end
% if ok
%     D1.channels.name=Name;
% end
        

%Plot glass brain
switch D1.Atlas
    case 'Human'
        P1=spm('dir');
    case {'Rat','Mouse','Macaque'}
        P1=fullfile(spm('dir'),'atlas5',lower(D1.Atlas));
end
path(P1,path);
% global defaults
% spm_defaults;
spm_figure('Clear',Fgraph)
hMIPax = axes('Parent',Fgraph,'Position',[0.13 0.52 0.715 0.468],'Visible','off');
hMIPax = spm_mip_ui([],zeros(3,0),eye(4),[2 2 2]',hMIPax);
h.hMIPax=hMIPax;

%-Axis offsets for 3d MIPs:
%=======================================================================
try
    Species=D1.Atlas;
catch
    Species=spm_input('Species ','+1','Human|Rat|Mouse');
end
%-MIP pane dimensions and Talairach origin offsets
%-See spm_project.c for derivation
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
switch Species
    case{'Human'}
        %             DXYZ = [182 218 182];
        %             CXYZ = [091 127 073];
        DXYZ = [192 218 182];
        CXYZ = [096 132 078];
    case{'Rat','Mouse'}
        DXYZ = [189 211 175];
        CXYZ = [106 167 141];
%         Pos=10*Pos;
end
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;
h.Po=Po;

h.D1=D1;
h.offset1=1;
h.offset2=1;
h.pvalue=min(D1.ca.p)+eps;

if length(D1.ca.Time)>1
%     %Slider for time
%     handles.timeslider=uicontrol(Fgraph, 'Tag', 'timeslider', 'Style', 'slider',...
%         'Min', 1, 'Max', D1.Nsamples, 'Value', 1, 'Units',...
%         'normalized', 'Position', [0.425 0.48 0.2 0.03],...
%         'SliderStep', [1/(D1.Nsamples-1) unique(max([1/(D1.Nsamples-1) 0.2]))],...
%         'TooltipString', 'Navigate through time',...
%         'Callback', @timeslider_update,...
%         'Parent', Fgraph, 'Interruptible', 'off');
%     uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
%         'String', 'Time [s]',...
%         'Position',[0.425 0.45 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', 'w');
%     uicontrol(Fgraph, 'Style', 'text', 'Tag', 'timedisp',...
%         'Units', 'normalized',...
%         'String', num2str(D1.time(get(handles.timeslider, 'Value'))),...
%         'Position',[0.425 0.511 0.2 0.029],...
%         'HorizontalAlignment', 'center', 'FontSize', FS1,...
%         'BackgroundColor', spm('Colour'));
    
    %p value update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'p value ',...
        'Position',[0.2 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'pvalue',...
        'Units', 'normalized',...
        'String', num2str(h.pvalue),...
        'Position',[0.3 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @pvalue_update,...
        'BackgroundColor', spm('Colour'));
    
    %range update
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] ',...
        'Position',[0.58 0.494 0.1 0.028],...
        'HorizontalAlignment', 'right', 'FontSize', FS2,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'rangetime',...
        'Units', 'normalized',...
        'String', num2str(h.RangeTime),...
        'Position',[0.68 0.494 0.1 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @rangetime_update,...
        'BackgroundColor', spm('Colour'));
    
    %Display Time
    uicontrol(Fgraph, 'Style', 'text', 'Units', 'normalized',...
        'String', 'Time [s] = ',...
        'Position',[0.55 0.55 0.15 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS1,...
        'BackgroundColor', 'w');
    uicontrol(Fgraph, 'Style', 'edit', 'Tag', 'currenttime',...
        'Units', 'normalized',...
        'String', num2str(h.D1.ca.Time(h.offset1)),...
        'Position',[0.69 0.55 0.08 0.028],...
        'HorizontalAlignment', 'left', 'FontSize', FS2,...
        'Callback', @currenttime_update,...
        'BackgroundColor', spm('Colour'));

    %Movie
    uicontrol(Fgraph, 'Style', 'pushbutton', 'Tag', 'movie',...
        'Units', 'normalized',...
        'String', 'Movie',...
        'Position',[0.9 0.55 0.08 0.028],...
        'HorizontalAlignment', 'center', 'FontSize', FS2,...
        'Callback', @makemovie,...
        'BackgroundColor', spm('Colour'));
    

    %Axes for ar data
    nchannels(D1)=length(D1.ca.channels);
    SizeAxeX=0.55/nchannels(D1);
    SizeAxeY=0.43/nchannels(D1);
    for i1=1:nchannels(D1)
        for i2=1:nchannels(D1)
            h.Heegaxes1(i1,i2) = axes('Position',[0.07+(i1-1)*SizeAxeX 0.02+(nchannels(D1)-i2)*SizeAxeY 0.95*SizeAxeX 0.95*SizeAxeY],...
                'Parent', Fgraph);
%             axes(h.Heegaxes1(i1,i2));
%             area(D1.ca.Time,squeeze(DataPlot1(i1,i2,:)));
%             h.eventaxes1(i1,i2)=line([D1.time(h.offset1) D1.time(h.offset1)],[0 1],'color','m','linewidth',1,'Visible','off');
        end
    end
    h.SizeAxeX=SizeAxeX;
    h.SizeAxeY=SizeAxeY;

    
    %Axes for zoom
    h.Heegaxes2(1) = axes('Position',[0.66 0.34 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(1));
%     h.eventaxes2{1} = line([D1.time(h.offset2) D1.time(h.offset2)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(2) = axes('Position',[0.66 0.19 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(2));
%     h.eventaxes2{2} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    h.Heegaxes2(3) = axes('Position',[0.66 0.04 0.32 0.1],'Parent', Fgraph);
%     axes(h.Heegaxes2(3));
%     h.eventaxes2{3} = line([D1.time(h.offset1) D1.time(h.offset1)-10],[mis mas],'color','m','linewidth',1,'Visible','off');
    
% 
%         if i1<=nchannels(D1)
%             h.Heegplot1{i1}=plot(D1.time,D1.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D1.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot1{i1}=plot(D1.time,D1.data(:,:));
%             hold on
%             h.Heegplot1{i1}=[h.Heegplot1{i1};plot(D1.time,mean(D1.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Coupling strength','FontSize', FS1)
% %         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes1(i1),'Visible','off')
%         set(h.Heegplot1{i1},'Visible','off')
%         h.eventaxes1{i1,1}=line([D1.time(1) D1.time(end)],[Zmax Zmax],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,2}=line([D1.time(1) D1.time(end)],[Zmin Zmin],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes1{i1,3}=line([D1.time(h.offset) D1.time(h.offset)],[mis mas],'color','k','linewidth',1,'Visible','off');
%     end

%     for i1=1:D2.Nchannels+1
%         h.Heegaxes2(i1) = axes('Position',...
%             [0.13 0.05 0.83 0.15],...
%             'Parent', Fgraph);
%         axes(h.Heegaxes2(i1));
%         if i1<=D2.Nchannels
%             h.Heegplot2{i1}=plot(D2.time,D2.data(i1,:),'k');
%             ylabel([sprintf('%d: ',i1) D2.channels.name{i1}], 'Interpreter', 'tex', 'FontSize', FS2);
%         else
%             h.Heegplot2{i1}=plot(D2.time,D2.data(:,:));
%             hold on
%             h.Heegplot2{i1}=[h.Heegplot2{i1};plot(D2.time,mean(D2.data(:,:)),'LineWidth',4,'Color','k')];
%             hold off
%             ylabel('All channels', 'Interpreter', 'tex', 'FontSize', FS2);
%         end
%         axis tight
%         title('Delay','FontSize', FS1)
%         xlabel('Time [s]', 'Interpreter', 'tex', 'FontSize', FS1);
%         set(h.Heegaxes2(i1),'Visible','off')
%         set(h.Heegplot2{i1},'Visible','off')
%         h.eventaxes2{i1,1}=line([D2.time(1) D2.time(end)],[Zd Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,2}=line([D2.time(1) D2.time(end)],[-Zd -Zd],'color','k','linewidth',1,'linestyle',':','Visible','off');
%         h.eventaxes2{i1,3}=line([D2.time(h.offset) D2.time(h.offset)],[-md md],'color','k','linewidth',1,'Visible','off');
%     end
                
    h.PairSelect=[1 1];
    guidata(Fgraph,h);
    [h.Threshold1,h.Threshold2]=pvalue_compute;
    guidata(Fgraph,h);
    h=draw_axes_gs;    
        
end


colordot='k';
axes(hMIPax)
for i1=1:size(Pos,1)
    hold on
    h1=plot(Po(1) + Pos(i1,2),Po(2) + Pos(i1,1),'.');
    h2=plot(Po(1) + Pos(i1,2),Po(3) - Pos(i1,3),'.');
    h3=plot(Po(4) + Pos(i1,1),Po(3) - Pos(i1,3),'.');
    hold off
    set(h1,'MarkerSize',15,'color',colordot);
    set(h2,'MarkerSize',15,'color',colordot);
    set(h3,'MarkerSize',15,'color',colordot);
end
% hquiv=zeros(6,size(D1.ca.CM,2));
% hold on
% for i1=1:size(hquiv,2)
%     hquiv(1,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(2,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(3,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Visible','off','Parent',Fgraph);
%     hquiv(4,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(2) + Pos(D1.ca.CM(2,i1),1),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(5,i1)=quiver(Po(1) + Pos(D1.ca.CM(1,i1),2),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),2)-Pos(D1.ca.CM(1,i1),2),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
%     hquiv(6,i1)=quiver(Po(4) + Pos(D1.ca.CM(1,i1),1),Po(3) - Pos(D1.ca.CM(2,i1),3),Pos(D1.ca.CM(2,i1),1)-Pos(D1.ca.CM(1,i1),1),Pos(D1.ca.CM(2,i1),3)-Pos(D1.ca.CM(1,i1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Visible','off','Parent',Fgraph);
% end
% hold off

% h.hquiv=hquiv;
% store handles
guidata(Fgraph, h);

% draw_lines;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timeslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

%save
h.offset=ind;

% make plots
for i1=1:size(h.eventaxes1,1)
    set(h.eventaxes1{i1,3},'XData',[h.D1.time(h.offset) h.D1.time(h.offset)]);
    set(h.eventaxes2{i1,3},'XData',[h.D1.time(h.offset) h.D1.time(h.offset)]);
end
guidata(F,h);
draw_lines;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delayslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

%save
h.Zd=1e-3*ind;

%disp value
set(handles.delaydisp,'String',num2str(h.Zd))

% make plots
for i1=1:size(h.eventaxes2,1)
    set(h.eventaxes2{i1,1},'YData',[h.Zd h.Zd]);
    set(h.eventaxes2{i1,2},'YData',[-h.Zd -h.Zd]);
end
guidata(F,h);
draw_lines

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function couplingmaxslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

%save
h.Zmax=1e-3*ind;

%disp value
set(handles.couplingmaxdisp,'String',num2str(h.Zmax))

% make plots
for i1=1:size(h.eventaxes1,1)
    set(h.eventaxes1{i1,1},'YData',[h.Zmax h.Zmax]);
end
guidata(F,h);
draw_lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function couplingminslider_update(hObject, events)
% update called from zoom buttons

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

% slider value
ind = round(get(hObject, 'Value'));
set(hObject, 'Value', ind);

%save
h.Zmin=1e-3*ind;

%disp value
set(handles.couplingmindisp,'String',num2str(h.Zmin))

% make plots
for i1=1:size(h.eventaxes1,1)
    set(h.eventaxes1{i1,2},'YData',[h.Zmin h.Zmin]);
end
guidata(F,h);
draw_lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_axes(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
% handles = guihandles(F);
for i1=1:length(h.Heegaxes1)
    set(h.Heegaxes1(i1),'Visible','off')
    set(h.Heegaxes2(i1),'Visible','off')
    set(h.Heegplot1{i1},'Visible','off')
    set(h.Heegplot2{i1},'Visible','off')
    set(h.eventaxes1{i1,1},'Visible','off')
    set(h.eventaxes1{i1,2},'Visible','off')
    set(h.eventaxes1{i1,3},'Visible','off')
    set(h.eventaxes2{i1,1},'Visible','off')
    set(h.eventaxes2{i1,2},'Visible','off')
    set(h.eventaxes2{i1,3},'Visible','off')
end
tmp=get(h.plotelec,'Value');
set(h.Heegaxes1(tmp),'Visible','on')
set(h.Heegaxes2(tmp),'Visible','on')
set(h.Heegplot1{tmp},'Visible','on')
set(h.Heegplot2{tmp},'Visible','on')
set(h.eventaxes1{tmp,1},'Visible','on')
set(h.eventaxes1{tmp,2},'Visible','on')
set(h.eventaxes1{tmp,3},'Visible','on')
set(h.eventaxes2{tmp,1},'Visible','on')
set(h.eventaxes2{tmp,2},'Visible','on')
set(h.eventaxes2{tmp,3},'Visible','on')
guidata(F,h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_data(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
tmp=get(h.load_data,'Value');
switch tmp
    case 1
        h.DataType='PSD';
        DataPlot1=h.D1.ca.S;
    case 2
        h.DataType='DTF';
        DataPlot1=h.D1.ca.DTF;
    case 3
        h.DataType='ffDTF';
        DataPlot1=h.D1.ca.ffDTF;
    case 4
        h.DataType='dDTF';
        DataPlot1=h.D1.ca.dDTF;
    case 5
        h.DataType='GdDTF';
        DataPlot1=h.D1.ca.GdDTF;
    case 6
        h.DataType='PDC';
        DataPlot1=h.D1.ca.PDC;
    case 7
        h.DataType='ffPDC';
        DataPlot1=h.D1.ca.ffPDC;
    case 8
        h.DataType='GPDC';
        DataPlot1=h.D1.ca.GPDC;
    case 9
        h.DataType='GffPDC';
        DataPlot1=h.D1.ca.GffPDC;
end
switch h.Mode
    case 'Time'
        DataPlot1=squeeze(mean(DataPlot1(:,:,:,h.RangeFreqIndex),4));
    case 'Freq'
        DataPlot1=squeeze(mean(DataPlot1(:,:,h.RangeTimeIndex,:),3));
end
DataPlot1=permute(DataPlot1,[2 1 3]);
DataPlot2=DataPlot1-permute(DataPlot1,[2 1 3]);       
h.DataPlot1=DataPlot1;
h.DataPlot2=DataPlot2;
guidata(F,h);
[h.Threshold1,h.Threshold2]=pvalue_compute;
guidata(F,h);
draw_axes_ar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=draw_axes_ar(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
SizeAxeX=h.SizeAxeX;
SizeAxeY=h.SizeAxeY;
switch h.Mode
    case 'Time'
        XAxes=h.D1.ca.Time;
        h.RangeIndex=h.RangeTimeIndex;
    case 'Freq'
        XAxes=h.D1.ca.Freq;
        h.RangeIndex=h.RangeFreqIndex;
end        
for i1=1:size(h.Heegaxes1,1)
    for i2=1:size(h.Heegaxes1,2)
        axes(h.Heegaxes1(i1,i2));
        if h.PairSelect(1)==i1&h.PairSelect(2)==i2
            tmp=area(XAxes(h.RangeIndex),ones(1,length(h.RangeIndex)));
            set(tmp,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
            hold on
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
            hold off
        elseif h.PairSelect(1)==i2&h.PairSelect(2)==i1
            tmp=area(XAxes(h.RangeIndex),ones(1,length(h.RangeIndex)));
            set(tmp,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9])
            hold on
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
            hold off
        else
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
        end
%         hold on
%         line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold1(i1,i2) h.Threshold1(i1,i2)],'Color','m')
%         hold off
        hold on
        tmp=squeeze(h.DataPlot1(i1,i2,h.RangeIndex));
        tmp(find(tmp<h.Threshold1(i1,i2)))=0;
        tmp=area(XAxes(h.RangeIndex),tmp);
        set(tmp,'FaceColor','r','EdgeColor','r')
        line([min(XAxes) max(XAxes)],[0 0],'Color','k')
        line([XAxes(h.offset1) XAxes(h.offset1)],[0 1],'color','m','linewidth',1,'Tag','Time');
        hold off
        axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) 0 max(h.DataPlot1(:))])
%         if h.PairSelect(1)==i1&h.PairSelect(2)==i2
%             set(tmp,'FaceColor','b','EdgeColor','b')
%         end
        set(h.Heegaxes1(i1,i2),'ButtonDownFcn','ImaGIN_DispCausality(''draw_axes_ar_button'')','Tag',num2str([i1 i2]));
        if i1>1
            set(h.Heegaxes1(i1,i2),'YTickLabel','')
        end
        if i2<size(h.Heegaxes1,2)
            set(h.Heegaxes1(i1,i2),'XTickLabel','')
        end
        if i2==1
            title(h.D1.chanlabels{i1})
        end
        if i1==1
            ylabel(h.D1.chanlabels{i2})
        end
    end
end
axes(h.Heegaxes2(1));
tmp=squeeze(h.D1(unique(h.PairSelect),:));
plot(h.D1.time,tmp);
title('Original EEG','FontSize',spm('FontSize', 8))
if strcmp(h.Mode,'Time')
    axis([min(h.D1.ca.Time(h.RangeIndex)) max(h.D1.ca.Time(h.RangeIndex)) min(tmp(:)) max(tmp(:))])
else
    axis([min(h.D1.ca.Time) max(h.D1.ca.Time) min(tmp(:)) max(tmp(:))])
end
hold on
line([h.D1.time(h.offset2) h.D1.time(h.offset2)],[min(tmp(:)) max(tmp(:))],'color','m','linewidth',1,'Tag','Time');
hold off
tmp=legend(h.D1.chanlabels{unique(h.PairSelect)});
legend('boxoff')
set(tmp,'FontSize',spm('FontSize',6),'Location','SouthWest')
set(h.Heegaxes2(1),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')
axes(h.Heegaxes2(2));
tmpplot=squeeze(h.DataPlot1(h.PairSelect(1),h.PairSelect(2),h.RangeIndex));
tmp=area(XAxes(h.RangeIndex),tmpplot);
% set(tmp,'FaceColor','b','EdgeColor','b')
hold on
tmpplot(find(tmpplot<h.Threshold1(h.PairSelect(1),h.PairSelect(2))))=0;
tmp=area(XAxes(h.RangeIndex),tmpplot);
set(tmp,'FaceColor','r','EdgeColor','r')
line([min(XAxes) max(XAxes)],[0 0],'Color','k')
line([XAxes(h.offset1) XAxes(h.offset1)],[0 1],'color','m','linewidth',1,'Tag','Time');
hold off
% hold on
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold1(h.PairSelect(1),h.PairSelect(2)) h.Threshold1(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% hold off
title([h.DataType ' from ' h.D1.chanlabels{unique(h.PairSelect(1))} ' to ' h.D1.chanlabels{unique(h.PairSelect(2))}],'FontSize',spm('FontSize', 8))
axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) 0 max(h.DataPlot1(:))])
set(h.Heegaxes2(2),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')
axes(h.Heegaxes2(3));
tmpplot=squeeze(h.DataPlot2(h.PairSelect(1),h.PairSelect(2),h.RangeIndex));
area(XAxes(h.RangeIndex),tmpplot)
hold on
tmp=tmpplot;
tmp(find(abs(tmp)<h.Threshold2(h.PairSelect(1),h.PairSelect(2))))=0;
tmp=area(XAxes(h.RangeIndex),tmp);
set(tmp,'FaceColor','r','EdgeColor','r')
line([min(XAxes) max(XAxes)],[0 0],'Color','k')
line([XAxes(h.offset1) XAxes(h.offset1)],[-1 1],'color','m','linewidth',1,'Tag','Time');
hold off
% hold on
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold2(h.PairSelect(1),h.PairSelect(2)) h.Threshold2(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],-[h.Threshold2(h.PairSelect(1),h.PairSelect(2)) h.Threshold2(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% hold off
title(['Difference of ' h.DataType ' from ' h.D1.chanlabels{unique(h.PairSelect(1))} ' to ' h.D1.chanlabels{unique(h.PairSelect(2))}],'FontSize',spm('FontSize', 8))
xlabel('Time [s]','FontSize',spm('FontSize', 8))
if sum(abs(tmpplot))>0
    axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) -max(abs(tmpplot)) max(abs(tmpplot))])
else
    xlim([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex))])
end
set(h.Heegaxes2(3),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')

guidata(F,h);

h=draw_lines_ar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=draw_axes_gs(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
SizeAxeX=h.SizeAxeX;
SizeAxeY=h.SizeAxeY;
XAxes=h.D1.ca.Time;
h.RangeIndex=h.RangeTimeIndex;
for i1=1:size(h.Heegaxes1,1)
    for i2=1:size(h.Heegaxes1,2)
        axes(h.Heegaxes1(i1,i2));
        if h.PairSelect(1)==i1&h.PairSelect(2)==i2
            tmp=area(XAxes(h.RangeIndex),ones(1,length(h.RangeIndex)));
            set(tmp,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
            hold on
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
            hold off
        elseif h.PairSelect(1)==i2&h.PairSelect(2)==i1
            tmp=area(XAxes(h.RangeIndex),ones(1,length(h.RangeIndex)));
            set(tmp,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9])
            hold on
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
            hold off
        else
            tmp=area(XAxes(h.RangeIndex),squeeze(h.DataPlot1(i1,i2,h.RangeIndex)));
        end
%         hold on
%         line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold1(i1,i2) h.Threshold1(i1,i2)],'Color','m')
%         hold off
        hold on
        tmp=squeeze(h.DataPlot1(i1,i2,h.RangeIndex));
        tmp(find(tmp<h.Threshold1(i1,i2)))=0;
        tmp=area(XAxes(h.RangeIndex),tmp);
        set(tmp,'FaceColor','r','EdgeColor','r')
        line([min(XAxes) max(XAxes)],[0 0],'Color','k')
        line([XAxes(h.offset1) XAxes(h.offset1)],[0 1],'color','m','linewidth',1,'Tag','Time');
        hold off
        axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) min(h.DataPlot1(:)) max(h.DataPlot1(:))])
%         if h.PairSelect(1)==i1&h.PairSelect(2)==i2
%             set(tmp,'FaceColor','b','EdgeColor','b')
%         end
        set(h.Heegaxes1(i1,i2),'ButtonDownFcn','ImaGIN_DispCausality(''draw_axes_ar_button'')','Tag',num2str([i1 i2]));
        if i1>1
            set(h.Heegaxes1(i1,i2),'YTickLabel','')
        end
        if i2<size(h.Heegaxes1,2)
            set(h.Heegaxes1(i1,i2),'XTickLabel','')
        end
        if i2==1
            title(h.D1.chanlabels{i1})
        end
        if i1==1
            ylabel(h.D1.chanlabels{i2})
        end
    end
end
% axes(h.Heegaxes2(1));
% tmp=squeeze(h.D1.data(unique(h.PairSelect),:));
% plot(h.D1.time,tmp);
% title('Original EEG','FontSize',spm('FontSize', 8))
% axis([min(h.D1.ca.Time(h.RangeIndex)) max(h.D1.ca.Time(h.RangeIndex)) min(tmp(:)) max(tmp(:))])
% hold on
% line([h.D1.time(h.offset2) h.D1.time(h.offset2)],[min(tmp(:)) max(tmp(:))],'color','m','linewidth',1,'Tag','Time');
% hold off
% tmp=legend(h.D1.chanlabels{unique(h.PairSelect)});
% legend('boxoff')
% set(tmp,'FontSize',spm('FontSize',6),'Location','SouthWest')
% set(h.Heegaxes2(1),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')
axes(h.Heegaxes2(2));
tmpplot=squeeze(h.DataPlot1(h.PairSelect(1),h.PairSelect(2),h.RangeIndex));
tmp=area(XAxes(h.RangeIndex),tmpplot);
% set(tmp,'FaceColor','b','EdgeColor','b')
hold on
tmpplot(find(tmpplot<h.Threshold1(h.PairSelect(1),h.PairSelect(2))))=0;
tmp=area(XAxes(h.RangeIndex),tmpplot);
set(tmp,'FaceColor','r','EdgeColor','r')
line([min(XAxes) max(XAxes)],[0 0],'Color','k')
line([XAxes(h.offset1) XAxes(h.offset1)],[0 1],'color','m','linewidth',1,'Tag','Time');
hold off
% hold on
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold1(h.PairSelect(1),h.PairSelect(2)) h.Threshold1(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% hold off
title(['GS from ' h.D1.chanlabels{unique(h.PairSelect(1))} ' to ' h.D1.chanlabels{unique(h.PairSelect(2))}],'FontSize',spm('FontSize', 8))
axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) min(h.DataPlot1(:)) max(h.DataPlot1(:))])
set(h.Heegaxes2(2),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')
axes(h.Heegaxes2(3));
tmpplot=squeeze(h.DataPlot2(h.PairSelect(1),h.PairSelect(2),h.RangeIndex));
area(XAxes(h.RangeIndex),tmpplot)
hold on
tmp=tmpplot;
tmp(find(abs(tmp)<h.Threshold2(h.PairSelect(1),h.PairSelect(2))))=0;
tmp=area(XAxes(h.RangeIndex),tmp);
set(tmp,'FaceColor','r','EdgeColor','r')
line([min(XAxes) max(XAxes)],[0 0],'Color','k')
line([XAxes(h.offset1) XAxes(h.offset1)],[-1 1],'color','m','linewidth',1,'Tag','Time');
hold off
% hold on
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],[h.Threshold2(h.PairSelect(1),h.PairSelect(2)) h.Threshold2(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% line([min(h.D1.ca.Time) max(h.D1.ca.Time)],-[h.Threshold2(h.PairSelect(1),h.PairSelect(2)) h.Threshold2(h.PairSelect(1),h.PairSelect(2))],'Color','m')
% hold off
title(['Difference of GS from ' h.D1.chanlabels{unique(h.PairSelect(1))} ' to ' h.D1.chanlabels{unique(h.PairSelect(2))}],'FontSize',spm('FontSize', 8))
xlabel('Time [s]','FontSize',spm('FontSize', 8))
if sum(abs(tmpplot))>0
    axis([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex)) -max(abs(tmpplot)) max(abs(tmpplot))])
else
    xlim([min(XAxes(h.RangeIndex)) max(XAxes(h.RangeIndex))])
end
set(h.Heegaxes2(3),'ButtonDownFcn','ImaGIN_DispCausality(''draw_time_button'')')

guidata(F,h);

h=draw_lines_ar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_axes_ar_button(varargin)

F=spm_figure('GetWin','Graphics');
h=guidata(F);
h.PairSelect=str2num(get(gco,'Tag'));
guidata(F,h);
if isfield(h.D1.ca,'GS')||isfield(h.D1.ca,'H2')
    draw_axes_gs;
else
    draw_axes_ar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_time_button(varargin)

F=spm_figure('GetWin','Graphics');
h=guidata(F);
[X,Y]=ginput(1);

switch h.Mode
    case 'Time'
        XAxes=h.D1.ca.Time;
        h.offset2=find(abs(h.D1.time-X)==min(abs(h.D1.time-X)));
    case 'Freq'
        XAxes=h.D1.ca.Freq;
        h.offset2=1;
end        
h.offset1=find(abs(XAxes-X)==min(abs(XAxes-X)));

for i1=1:size(h.Heegaxes1,1)
    for i2=1:size(h.Heegaxes1,2)
        axes(h.Heegaxes1(i1,i2));
        tmp=get(gca,'Children');
        for i3=1:length(tmp)
            tmp2=get(tmp(i3),'Tag');
            if strcmp(tmp2,'Time')
                set(tmp(i3),'XData',[XAxes(h.offset1) XAxes(h.offset1)])
                break
            end
        end
    end
end
for i1=1:3
    axes(h.Heegaxes2(i1));
    tmp=get(gca,'Children');
    for i3=1:length(tmp)
        tmp2=get(tmp(i3),'Tag');
        if strcmp(tmp2,'Time')
            if i1==1
                switch h.Mode
                    case 'Time'
                        set(tmp(i3),'XData',[h.D1.time(h.offset2) h.D1.time(h.offset2)])
                    case 'Freq'
                        set(tmp(i3),'XData',[h.D1.time(h.offset2)-10 h.D1.time(h.offset2)-10])
                end
                tmp=legend(h.D1.chanlabels{unique(h.PairSelect)});
                legend('boxoff')
                set(tmp,'FontSize',spm('FontSize',6),'Location','SouthWest')
            else
                set(tmp(i3),'XData',[XAxes(h.offset1) XAxes(h.offset1)])
            end
            break
        end
    end
end
handles=guihandles(F);
set(handles.currenttime,'String',num2str(h.D1.ca.Time(h.offset1)))
guidata(F,h);
draw_lines_ar;
% draw_axes_ar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function currenttime_update(varargin)

F  = spm_figure('GetWin','Graphics');

handles=guihandles(F);
t=str2num(get(handles.currenttime,'String'));

h = guidata(F);
if strcmp(h.Mode,'Time')
    XAxes=h.D1.ca.Time;
    if t<h.D1.ca.Time(1)
        t=h.D1.ca.Time(1);
        set(gco,'String',num2str(t))
    elseif t>h.D1.ca.Time(end)
        t=h.D1.ca.Time(end);
        set(gco,'String',num2str(t))
    end
    h.offset1=find(abs(h.D1.ca.Time-t)==min(abs(h.D1.ca.Time-t)));
    h.offset2=find(abs(h.D1.time-t)==min(abs(h.D1.time-t)));

    % draw_axes_ar;

    for i1=1:size(h.Heegaxes1,1)
        for i2=1:size(h.Heegaxes1,2)
            axes(h.Heegaxes1(i1,i2));
            tmp=get(gca,'Children');
            for i3=1:length(tmp)
                tmp2=get(tmp(i3),'Tag');
                if strcmp(tmp2,'Time')
                    set(tmp(i3),'XData',[XAxes(h.offset1) XAxes(h.offset1)])
                    break
                end
            end
        end
    end
    for i1=1:3
        axes(h.Heegaxes2(i1));
        tmp=get(gca,'Children');
        for i3=1:length(tmp)
            tmp2=get(tmp(i3),'Tag');
            if strcmp(tmp2,'Time')
                if i1==1
                    switch h.Mode
                        case 'Time'
                            set(tmp(i3),'XData',[h.D1.time(h.offset2) h.D1.time(h.offset2)])
                        case 'Freq'
                            set(tmp(i3),'XData',[h.D1.time(h.offset2)-10 h.D1.time(h.offset2)-10])
                    end
                    tmp=legend(h.D1.chanlabels{unique(h.PairSelect)});
                    legend('boxoff')
                    set(tmp,'FontSize',spm('FontSize',6),'Location','SouthWest')
                else
                    set(tmp(i3),'XData',[XAxes(h.offset1) XAxes(h.offset1)])
                end
                break
            end
        end
    end
    guidata(F,h);
    draw_lines_ar;
    % draw_axes_ar;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makemovie(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles=guihandles(F);
t0=get(handles.currenttime,'String');

direc = spm_select(1,'dir','Select directory in which movie is saved ');
filename=spm_input('Name of movie ', 1, 's');
filename=fullfile(direc,filename);
timerange=spm_input('Time range [s] ', '+1', 'r',[h.D1.ca.Time(1) h.D1.ca.Time(end)],2);   %time range
t=find(h.D1.ca.Time>=timerange(1)&h.D1.ca.Time<=timerange(2));
fps=spm_input('Frames per second ', '+1', 'i',3,1);   %time range
%loop sur les frames
for i1=1:length(t)
    set(handles.currenttime,'String',num2str(h.D1.ca.Time(t(i1))));
    currenttime_update;
    M(i1)=getframe(F,[90 515 457 384]);
end
movie2avi(M,filename,'FPS',fps);
set(handles.currenttime,'String',t0);
currenttime_update;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Threshold1,Threshold2]=pvalue_compute(varargin)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

if isfield(h.D1.ca,'GS_S')
    data=h.D1.ca.GS_S;
    data=permute(data,[2 1 3]);
    datadiff=abs(data-permute(data,[2 1 3]));
    data=sort(data,3,'descend');
    datadiff=sort(datadiff,3,'descend');
    tmp=max(find(h.D1.ca.p<=h.pvalue));
    if isempty(tmp)
        tmp=1;
    end
    Threshold1=squeeze(data(:,:,tmp));
    Threshold2=squeeze(datadiff(:,:,tmp));
elseif isfield(h.D1.ca,'H2_S')
    data=h.D1.ca.H2_S;
    data=permute(data,[2 1 3]);
    datadiff=abs(data-permute(data,[2 1 3]));
    data=sort(data,3,'descend');
    datadiff=sort(datadiff,3,'descend');
    tmp=max(find(h.D1.ca.p<=h.pvalue));
    if isempty(tmp)
        tmp=1;
    end
    Threshold1=squeeze(data(:,:,tmp));
    Threshold2=squeeze(datadiff(:,:,tmp));
else


    switch h.DataType
        case 'PSD'
            data=h.D1.ca.S_S;
        case 'DTF'
            data=h.D1.ca.DTF_S;
        case 'ffDTF'
            data=h.D1.ca.ffDTF_S;
        case 'dDTF'
            data=h.D1.ca.dDTF_S;
        case 'GdDTF'
            data=h.D1.ca.GdDTF_S;
        case 'PDC'
            data=h.D1.ca.PDC_S;
        case 'ffPDC'
            data=h.D1.ca.ffPDC_S;
        case 'GPDC'
            data=h.D1.ca.GPDC_S;
        case 'GffPDC'
            data=h.D1.ca.GffPDC_S;
    end
    switch h.Mode
        case 'Time'
            data=squeeze(mean(data(:,:,:,h.RangeFreqIndex),4));
            data=permute(data,[2 1 3]);
            datadiff=abs(data-permute(data,[2 1 3]));
            data=sort(data,3,'descend');
            datadiff=sort(datadiff,3,'descend');
            tmp=max(find(h.D1.ca.p<=h.pvalue));
            if isempty(tmp)
                tmp=1;
            end
            Threshold1=squeeze(data(:,:,tmp));
            Threshold2=squeeze(datadiff(:,:,tmp));
        case 'Freq'
            data=sort(data,3,'descend');
            datadiff=abs(data-permute(data,[2 1 3 4]));
            datadiff=sort(datadiff,3,'descend');
            data=squeeze(mean(data(:,:,:,:),3));
            tmp=max(find(h.D1.ca.p<=h.pvalue));
            if isempty(tmp)
                tmp=1;
            end
            %         Threshold1=squeeze(data(:,:,tmp,:));
            %         Threshold2=squeeze(datadiff(:,:,tmp,:));
            Threshold1=max(max(data,[],4),[],3);    %No significance computed
            Threshold2=max(max(datadiff,[],4),[],3); %No significance computed
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pvalue_update(varargin)

p=str2num(get(gco,'String'));
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
if p<h.D1.ca.pmax
    p=h.D1.ca.pmax;
    set(gco,'String',num2str(p))
end
h.pvalue=p;
guidata(F,h);
[h.Threshold1,h.Threshold2]=pvalue_compute;
guidata(F,h);
if isfield(h.D1.ca,'GS_S')||isfield(h.D1.ca,'H2_S')
    draw_axes_gs;
else
    draw_axes_ar;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rangetime_update(varargin)

Range=str2num(get(gco,'String'));
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
XAxes=h.D1.ca.Time;
if min(Range)<min(XAxes)
    Range(1)=min(XAxes);
    set(gco,'String',num2str(Range))
end
if max(Range)>max(XAxes)
    Range(2)=max(XAxes);
    set(gco,'String',num2str(Range))
end
h.RangeTime=Range;
h.RangeTimeIndex=[find(abs(XAxes-Range(1))==min(abs(XAxes-Range(1)))):find(abs(XAxes-Range(2))==min(abs(XAxes-Range(2))))];
guidata(F,h);
switch h.Mode
    case 'Time'
        if isfield(h.D1.ca,'GS')||isfield(h.D1.ca,'H2')
            draw_axes_gs;
        else
            draw_axes_ar;
        end
    case 'Freq'
        load_data;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rangefreq_update(varargin)

Range=str2num(get(gco,'String'));
F  = spm_figure('GetWin','Graphics');
h = guidata(F);
XAxes=h.D1.ca.Freq;
if min(Range)<min(XAxes)
    Range(1)=min(XAxes);
    set(gco,'String',num2str(Range))
end
if max(Range)>max(XAxes)
    Range(2)=max(XAxes);
    set(gco,'String',num2str(Range))
end
h.RangeFreq=Range;
h.RangeFreqIndex=[find(abs(XAxes-Range(1))==min(abs(XAxes-Range(1)))):find(abs(XAxes-Range(2))==min(abs(XAxes-Range(2))))];
guidata(F,h);
switch h.Mode
    case 'Time'
        load_data;
    case 'Freq'
        draw_axes_ar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function change_mode(varargin)

F  = spm_figure('GetWin','Graphics');
handles=guihandles(F);
h = guidata(F);
tmp=get(h.change_mode,'Value');
switch tmp
    case 1
        h.Mode='Time';
        X=h.RangeTime(1);
        XAxes=h.D1.ca.Time;
        h.offset2=find(abs(h.D1.time-X)==min(abs(h.D1.time-X)));
    case 2
        h.Mode='Freq';
        X=h.RangeFreq(1);
        XAxes=h.D1.ca.Freq;
        h.offset2=1;
end
h.offset1=find(abs(XAxes-X)==min(abs(XAxes-X)));
% set(handles.range,'String',num2str(h.Range))
% DataPlot2=DataPlot1-permute(DataPlot1,[2 1 3]);       
% h.DataPlot1=DataPlot1;
% h.DataPlot2=DataPlot2;
guidata(F,h);
load_data;
% [h.Threshold1,h.Threshold2]=pvalue_compute;
% guidata(F,h);
% draw_axes_ar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_lines

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

D1=h.D1;
D2=h.D2;
Zmax=h.Zmax;
Zmin=h.Zmin;
Zd=h.Zd;
hMIPax=handles.hMIPax;
Po=h.Po;
Pos=h.Pos;
offset=h.offset;

PairPos1=find(D1(:,offset)>Zmax);
PairNeg1=find(D1(:,offset)<Zmin);
% PairPos2=find(D2.data(:)>0);
% PairNeg2=find(D2.data(:)<0);
M=D1.ca.CM;

% colordot='k';
colorline1='r';
colorline2='b';
% axes(hMIPax)
% for i1=1:size(Pos,1)
%     hold on
%     h1=plot(Po(1) + Pos(i1,2),Po(2) + Pos(i1,1),'.');
%     h2=plot(Po(1) + Pos(i1,2),Po(3) - Pos(i1,3),'.');
%     h3=plot(Po(4) + Pos(i1,1),Po(3) - Pos(i1,3),'.');
%     hold off
%     set(h1,'MarkerSize',15,'color',colordot);
%     set(h2,'MarkerSize',15,'color',colordot);
%     set(h3,'MarkerSize',15,'color',colordot);
% end
% set(hquiv,'Visible','on')
try
    hquiv=h.hquiv;
%     set(hquiv(find(hquiv)),'Visible','off')
    delete(hquiv(find(hquiv)))
end
hquiv=zeros(6,size(D1.ca.CM,2));
for i1=1:length(PairPos1)
    ind=M(:,PairPos1(i1));
    if D2.data(PairPos1(i1),offset)>0
        ind=flipdim(ind,1);
    end
%     scribe.arrow('X',Po(1) + Pos(ind,2),'Y',Po(2) + Pos(ind,1),'color',colorline1,'Parent',hMIPax,'HeadStyle','vback2')
%         set(h,'HeadLength',100,'HeadWidth',100)
    hold on
    hquiv(1,i1)=quiver(Po(1) + Pos(ind(1),2),Po(2) + Pos(ind(1),1),Pos(ind(2),2)-Pos(ind(1),2),Pos(ind(2),1)-Pos(ind(1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
    if abs(D2.data(PairPos1(i1),offset))<=Zd
        set(hquiv(1,i1),'ShowArrowHead','off')
    end
    hquiv(2,i1)=quiver(Po(1) + Pos(ind(1),2),Po(3) - Pos(ind(1),3),Pos(ind(2),2)-Pos(ind(1),2),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
    if abs(D2.data(PairPos1(i1),offset))<=Zd
        set(hquiv(2,i1),'ShowArrowHead','off')
    end
    hquiv(3,i1)=quiver(Po(4) + Pos(ind(1),1),Po(3) - Pos(ind(1),3),Pos(ind(2),1)-Pos(ind(1),1),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
    if abs(D2.data(PairPos1(i1),offset))<=Zd
        set(hquiv(3,i1),'ShowArrowHead','off')
    end
end
for i1=1:length(PairNeg1)
    ind=M(:,PairNeg1(i1));
    if D2.data(PairNeg1(i1),offset)>0
        ind=flipdim(ind,1);
    end
    hquiv(4,i1)=quiver(Po(1) + Pos(ind(1),2),Po(2) + Pos(ind(1),1),Pos(ind(2),2)-Pos(ind(1),2),Pos(ind(2),1)-Pos(ind(1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Parent',F);
    if abs(D2.data(PairNeg1(i1),offset))<=Zd
        set(hquiv(4,i1),'ShowArrowHead','off')
    end
    hquiv(5,i1)=quiver(Po(1) + Pos(ind(1),2),Po(3) - Pos(ind(1),3),Pos(ind(2),2)-Pos(ind(1),2),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Parent',F);
    if abs(D2.data(PairNeg1(i1),offset))<=Zd
        set(hquiv(5,i1),'ShowArrowHead','off')
    end
    hquiv(6,i1)=quiver(Po(4) + Pos(ind(1),1),Po(3) - Pos(ind(1),3),Pos(ind(2),1)-Pos(ind(1),1),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline2,'Parent',F);
    if abs(D2.data(PairNeg1(i1),offset))<=Zd
        set(hquiv(6,i1),'ShowArrowHead','off')
    end
end
hold off
h.hquiv=hquiv;
guidata(F,h);

set(handles.timedisp, 'String',num2str(D1.time(offset)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=draw_lines_ar

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
handles = guihandles(F);

hMIPax=handles.hMIPax;
Po=h.Po;
Pos=h.Pos;
offset=h.offset1;

[Select1,Select2]=find(squeeze(h.DataPlot2(:,:,h.offset1))>h.Threshold2);

colorline1='r';
axes(hMIPax)
try
    hquiv=h.hquiv;
    delete(hquiv(find(hquiv)))
end
hquiv=zeros(size(h.DataPlot2,1),size(h.DataPlot2,2));
for i1=1:length(Select1)
    ind=[Select1(i1) Select2(i1)];
%     ind=flipdim(ind,1);
    hold on
%     hquiv(1,i1)=quiver(Po(1) + Pos(ind(1),2),Po(2) + Pos(ind(1),1),Pos(ind(2),2)-Pos(ind(1),2),Pos(ind(2),1)-Pos(ind(1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
%     hquiv(2,i1)=quiver(Po(1) + Pos(ind(1),2),Po(3) - Pos(ind(1),3),Pos(ind(2),2)-Pos(ind(1),2),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
%     hquiv(3,i1)=quiver(Po(4) + Pos(ind(1),1),Po(3) - Pos(ind(1),3),Pos(ind(2),1)-Pos(ind(1),1),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1,'Parent',F);
    hquiv(1,i1)=quiver(Po(1) + Pos(ind(1),2),Po(2) + Pos(ind(1),1),Pos(ind(2),2)-Pos(ind(1),2),Pos(ind(2),1)-Pos(ind(1),1),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1);
    hquiv(2,i1)=quiver(Po(1) + Pos(ind(1),2),Po(3) - Pos(ind(1),3),Pos(ind(2),2)-Pos(ind(1),2),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1);
    hquiv(3,i1)=quiver(Po(4) + Pos(ind(1),1),Po(3) - Pos(ind(1),3),Pos(ind(2),1)-Pos(ind(1),1),-Pos(ind(2),3)+Pos(ind(1),3),0,'AutoScale','off','MaxHeadSize',5,'color',colorline1);
end
hold off
h.hquiv=hquiv;
guidata(F,h);

