function ImaGIN(cmd)
% IMAGIN Starts the ImaGIN toolbox GUI (spm toolbox).

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% Copyright (c)2000-2017 Inserm
% =============================================================================-
%
% Authors: Olivier David, 2010-2017
%          Francois Tadel, 2017

%-Set ImaGIN path
%-----------------------------------------------------------------------
% Make sure the main ImaGIN folder is at the top of the path
ImaGINdir = fileparts(which(mfilename));
addpath(ImaGINdir, '-BEGIN'); 
% List of folders to add
toolboxDir = {'external','toolbox'}; % in reverse order of priority
for i = 1:length(toolboxDir)
    nextdir = fullfile(ImaGINdir,toolboxDir{i});
    % Check that directory exist
    if ~isdir(nextdir)
        error(['Directory "' toolboxDir{i} '" does not exist in the ImaGIN path.' 10 ...
               'Please re-install ImaGIN.']);
    end
    % Recursive search for subfolders in each main folder
    P = genpath(nextdir);
    % Add directory and subdirectories
    addpath(P, '-BEGIN');
end


%-Switch for specific commands
%-----------------------------------------------------------------------
if (nargin >= 1) && ~isempty(cmd)
    switch lower(cmd)
        case 'deploy'
            addpath(fullfile(ImaGINdir, 'deploy'));
            ImaGIN_deploy
            return;
    end
end


%-Initializes SPM
%-----------------------------------------------------------------------
% Configure SPM
spm('Defaults','EEG')
% Vis='on';
Modalities = {'PET','FMRI','EEG'};
%-Close any existing 'Menu' 'Tag'ged windows
delete(spm_figure('FindWin','Menu'))

%-Get size and scalings and create Menu window
%-----------------------------------------------------------------------
WS   = spm('WinScale');				%-Window scaling factors
FS   = spm('FontSizes');			%-Scaled font sizes
PF   = spm_platform('fonts');			%-Font names (for this platform)
Rect = spm('WinSize','Menu','raw').*WS;		%-Menu window rectangle

[SPMver,SPMc] = spm('Ver','',1);

Fmenu = figure('IntegerHandle','off',...
	'Name',sprintf('%s%s',spm('ver'),spm('GetUser',' (%s)')),...
	'NumberTitle','off',...
	'Tag','Menu',...
	'Position',Rect,...
	'Resize','off',...
	'Color',[1 1 1]*.8,...
	'UserData',struct('SPMver',SPMver,'SPMc',SPMc),...
	'MenuBar','none',...
	'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(10),...
	'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(12),...
	'DefaultUicontrolInterruptible','on',...
	'Renderer','painters',...
	'Visible','on', ...
    'CloseRequestFcn', 'spm(''Quit'')');

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 145 380 295].*WS)

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 340 360 90].*WS)

uicontrol(Fmenu,'Style','Text',...
  'String','Data definition',...
	'Position',[025 410 350 020].*WS,...
	'ForegroundColor','k',...
    'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 295 360 40].*WS, 'Tag', 'EEG')

uicontrol(Fmenu,'Style','Text',...
	'String','Analysis',...
	'Position',[025 315 350 020].*WS,...
	'ForegroundColor','k',...
    'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...    
	'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 200 360 90].*WS,...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','Text',...
    'String','Display',...
	'Position',[025 270 350 020].*WS,...
	'ForegroundColor','k',...
	'FontSize',FS(10),...
    'FontAngle','Italic',...
	'HorizontalAlignment','Left',...
	'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 155 360 40].*WS,...
    'Tag', 'EEG')
uicontrol(Fmenu,'Style','Text',...
	'String','Inference',...
	'Position',[025 175 350 020].*WS,...
	'ForegroundColor','k',...
	'FontAngle','Italic',...
	'FontSize',FS(10),...
	'HorizontalAlignment','Left',...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','Text',...
	'String','SPM Intracerebral EEG Toolbox',...
	'ToolTipString','modality & defaults set for EEG/MEG',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontName',PF.times,'FontAngle','Italic','FontWeight','Bold',...
	'HorizontalAlignment','center',...
	'Position',[020 122 360 020].*WS,...
	'Tag','EEG','Visible','on')

% Utilities frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame',...
	'BackgroundColor',spm('Colour'),...
	'Position',[010 010 380 112].*WS)

uicontrol(Fmenu,'Style','Text','String',SPMc,...
	'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	'FontName',PF.times,...
	'FontSize',FS(10),...
	'HorizontalAlignment','center',...
	'Position',[020 002 360 008].*WS)

%-Objects with Callbacks - main spm_*_ui.m routines
%=======================================================================

%-Data
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Convert',...
    'Position',[035 380 100 030].*WS,...
	'ToolTipString','Convert data in SPM EEG format',...
	'FontSize',FS(10),...
    'UserData','ImaGIN_spm_eeg_converteeg2mat',...
	'CallBack','ImaGIN_spm_eeg_converteeg2mat;',...
    'Tag', 'EEG',...
    'Visible','on')

uicontrol(Fmenu,'String','Montage',...
    'Position',[035 345 100 030].*WS,...
	'ToolTipString','Create montage',...
	'FontSize',FS(10),...
    'UserData','ImaGIN_Montage',...
	'CallBack','ImaGIN_Montage(1);',...
    'Tag', 'EEG',...
    'Visible','on');
	
uicontrol(Fmenu,'String','Events',...
    'Position',[150 380 100 030].*WS,...
	'ToolTipString','Events',...
	'FontSize',FS(10),...
    'UserData','ImaGIN_Events',...
    'CallBack','ImaGIN_Events;',...
    'Tag', 'EEG',...
    'Visible','on');

uicontrol(Fmenu,'String','Time zero',...
    'Position',[150 345 100 030].*WS,...
    'ToolTipString','Define time origin compared to a pre-defined event',...
    'FontSize',FS(10),...
    'UserData','ImaGIN_TimeZero',...
    'CallBack','ImaGIN_TimeZero;',...
    'Tag', 'EEG',...
    'Visible','on');

uicontrol(Fmenu,'String','Epoch',...
    'Position',[265 380 100 030].*WS,...
	'ToolTipString','Select a time window according to pre-defined events and overwrite a new (smaller) data set',...
	'FontSize',FS(10),...
    'UserData','ImaGIN_Crop',...
	'CallBack','ImaGIN_Crop;',...
    'Tag', 'EEG',...
    'Visible','on');

uicontrol(Fmenu,'Style','PopUp',...
    'String','Other...|BP filter|Notch|rereference|downsample|time-frequency|grand mean|merge|3Dimages|average TF',...
	'Position',[265 345 100 030].*WS,...
	'ToolTipString','Other preprocessing functions',...
	'FontSize',FS(10),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_BPFilter;','ImaGIN_NotchFilter;','spm_eeg_rereference;','ImaGIN_spm_eeg_downsample;','spm_eeg_tf;',...
    'spm_eeg_grandmean;','spm_eeg_merge;','spm_eeg_make3dimage;','spm_eeg_average_TF;'},...
    'Tag', 'EEG')

% Spatial EEG
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','PopUp',...
    'String','FFT|  Compute|  Normalise|  Average',...
	'Position',[25 290 50 030].*WS,...
	'ToolTipString','Time frequency analysis (power and synchrony)',...
	'FontSize',FS(10),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_FFT;','ImaGIN_NormaliseFFT;','ImaGIN_AverageFFT;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','PopUp',...
    'String','Time Freq|  Compute|  Normalise|  Rescale|  Average',...
	'Position',[85 290 65 030].*WS,...
	'ToolTipString','Time frequency analysis (power and synchrony)',...
	'FontSize',FS(10),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_spm_eeg_tf;','ImaGIN_NormaliseTF;','ImaGIN_RescaleTF;','ImaGIN_AverageTF;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','PopUp',...
    'String','Causality|  Compute|  Normalise|  Rescale|  Average',...
	'Position',[160 290 65 030].*WS,...
	'ToolTipString','Causality analysis',...
	'FontSize',FS(10),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_causality_compute;','ImaGIN_causality_normalise;','ImaGIN_causality_rescale;','ImaGIN_causality_average;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','PopUp',...
    'String','Firing rate|  Compute|  Normalise|  Average',...
	'Position',[235 290 65 030].*WS,...
	'ToolTipString','Firing rate',...
	'FontSize',FS(10),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_spike_compute;','ImaGIN_spike_normalise;','ImaGIN_spike_average;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','3D interpolation',...
    'Position',[310 300 65 030].*WS,...
	'ToolTipString', '3D interpolation of EEG intracerebral data',...
	'FontSize',FS(10),...
	'UserData','ImaGIN_spm_eeg_TF_images_3D;',...
	'CallBack','ImaGIN_spm_eeg_TF_images_3D;',...
	'Visible','on',...
    'Tag','EEG')

uicontrol(Fmenu,'String','2D interpolation',...
    'Position',[120 300 100 030].*WS,...
	'ToolTipString', '2D interpolation of EEG/MEG channel data',...
	'FontSize',FS(10),...
	'UserData','spm_eeg_TF_images;',...
	'CallBack','spm_eeg_TF_images;',...
	'Visible','off',...
    'Tag','EEG')

uicontrol(Fmenu,'String','Average TF',...
    'Position',[235 300 135 030].*WS,...
	'ToolTipString','3D source reconstruction of EEG/MEG data',...
	'FontSize',FS(10),...
    'UserData','',...
	'CallBack','spm_eeg_inv_imag_api;',...
	'Visible','off',...
    'Tag','EEG')

% Statistical EEG
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Montage',...
	'Position',[25 240 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','Display montage',...
	'UserData','ImaGIN_Montage',...
	'CallBack','ImaGIN_Montage(2);',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Data',...
	'Position',[25 205 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','Display data',...
	'UserData','ImaGIN_DispData',...
	'CallBack','ImaGIN_DispData;',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Time Frequency',...
	'Position',[145 240 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','Display time frequency analysis',...
	'UserData','ImaGIN_DispTF',...
	'CallBack','ImaGIN_DispTF',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Brain overlay',...
    'Position',[145 205 110 030].*WS,...
	'ToolTipString', '3D EEG intracerebral data overlaid on a brain',...
	'FontSize',FS(10),...
	'UserData','ImaGIN_Disp3D;',...
	'CallBack','ImaGIN_Disp3D;',...
	'Visible','on',...
    'Tag','EEG')

uicontrol(Fmenu,'String','Causality',...
	'Position',[265 240 110 030].*WS,...
	'FontSize',FS(10),...
	'ToolTipString','Display causality analysis',...
	'UserData','ImaGIN_DispCausality',...
	'CallBack','ImaGIN_DispCausality',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Results',...
	'Position',[130 160 130 030].*WS,...
	'ToolTipString','Inference and regional responses etc.',...
	'FontSize',FS(10),...
	'UserData','spm_results_ui',...
	'CallBack','[hReg,xSPM,SPM] = spm_results_ui;',...
    'Tag', 'EEG');

% -Utility buttons (first line)
% -----------------------------------------------------------------------
uicontrol(Fmenu,'Style','PopUp',...
    'String','Display...|images|M/EEG',...
	'Position',[020 088 082 024].*WS,...
	'ToolTipString','orthogonal sections and M/EEG display',...
	'FontSize',FS(9),...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{	'spm_jobman(''serial'','''',''jobs.util.disp'');','spm_eeg_review;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','Check Reg',...
	'Position',[112 088 083 024].*WS,...
	'ToolTipString','check image registration',...
	'FontSize',FS(9),...
	'UserData','jobs.util.checkreg',...
	'CallBack','spm_jobman(''serial'','''',''jobs.util.checkreg'');')

uicontrol(Fmenu,'Style','PopUp',...
	'String',Modalities,...
	'ToolTipString','change modality PET<->fMRI<->EEG',...
	'Tag','Modality',...
	'Position',[298 088 082 024].*WS,...
	'CallBack',[...
        'if get(gco,''Value'')==1;',...
        'spm ''PET'';',...
        'elseif get(gco,''Value'')==2;',...
        'spm ''FMRI'';',...
        'elseif get(gco,''Value'')==3;',...
        'spm ''EEG'';',...
        'end'],...
    'UserData','SPM_Modality',...
	'Interruptible','off')


%-Utility buttons (second line)
%-----------------------------------------------------------------------

uicontrol(Fmenu,'String','Seizure detect',...
	'Position',[020 054 083 024].*WS,...
	'ToolTipString','Detect seizure and write a time series in a new channel',...
	'FontSize',FS(9),...
	'UserData','ImaGIN_SeizureDetect',...
	'CallBack','ImaGIN_SeizureDetect;')

uicontrol(Fmenu,'String','Epileptogenicity',...
	'Position',[020 020 083 024].*WS,...
	'ToolTipString','Compute epileptogenicity on SEEG',...
	'FontSize',FS(9),...
	'UserData','ImaGIN_Epileptogenicity',...
	'CallBack','ImaGIN_Epileptogenicity;')

uicontrol(Fmenu,'String','Stim detect',...
	'Position',[112 054 083 024].*WS,...
	'ToolTipString','Detect stimulation and write events in a text file',...
	'FontSize',FS(9),...
	'UserData','ImaGIN_StimDetect',...
	'CallBack','ImaGIN_StimDetect;')

uicontrol(Fmenu,'String','Artefact correct',...
	'Position',[112 020 083 024].*WS,...
	'ToolTipString','Correct stimulation artefacts',...
	'FontSize',FS(9),...
	'UserData','ImaGIN_ArtefactCorrectionModel',...
	'CallBack','ImaGIN_ArtefactCorrectionModel;')

uicontrol(Fmenu,'String','Implantation',...
	'Position',[205 020 083 024].*WS,...
	'ToolTipString','Prepare montage from postop MRI',...
	'FontSize',FS(9),...
	'UserData','ImaGIN_MRImplantation',...
	'CallBack','ImaGIN_MRImplantation;')

% uicontrol(Fmenu,'String','Atlas connect',...
% 	'Position',[205 088 083 024].*WS,...
% 	'ToolTipString','Display DES connectivity atlas',...
% 	'FontSize',FS(9),...
% 	'UserData','ImaGIN_ConnectStim',...
% 	'CallBack','ImaGIN_ConnectStim;')

% uicontrol(Fmenu,'String','ImCalc',...
% 	'Position',[205 054 083 024].*WS,...
% 	'ToolTipString','image calculator',...
% 	'FontSize',FS(9),...
% 	'UserData','jobs.util.imcalc',...
% 	'CallBack','spm_jobman(''interactive'','''',''jobs.util.imcalc'');')

uicontrol(Fmenu,'String','Quit',...
	'Position',[298 020 082 024].*WS,...
	'ToolTipString','exit SPM',...
	'ForeGroundColor','r',...
	'Interruptible','off',...
	'UserData','QuitSPM',...
	'CallBack','spm(''Quit''), clear all')




