function ImaGIN(cmd)
% Starts the ImaGIN toolbox GUI (spm toolbox).

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
% Check installation of SPM
if ~exist('spm', 'file')
    disp([10 '====================================================================']);
    disp('ImaGIN requires SPM12:');
    disp(' - Download SPM: <a href="http://www.fil.ion.ucl.ac.uk/spm/software/download/">http://www.fil.ion.ucl.ac.uk/spm/software/download/</a>');
    disp(' - Add the SPM folder to your Matlab path');
    disp(['====================================================================' 10]);
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
	... 'DefaultTextFontName',PF.helvetica,...
	'DefaultTextFontSize',FS(8),...
	... 'DefaultUicontrolFontName',PF.helvetica,...
	'DefaultUicontrolFontSize',FS(8),...
	'DefaultUicontrolInterruptible','on',...
	'Renderer','painters',...
	'Visible','on', ...
    'CloseRequestFcn', 'spm(''Quit'')');

%-Frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame','BackgroundColor',spm('Colour'),...
	'Position',[010 112 380 328].*WS)

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 340 360 87].*WS)
uicontrol(Fmenu,'Style','Text',...
  'String','Data definition',...
	'Position',[025 406 350 020].*WS,...
	'ForegroundColor','k',...
    'FontAngle','Italic',...
	'HorizontalAlignment','Left',...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 254 360 82].*WS, 'Tag', 'EEG')
uicontrol(Fmenu,'Style','Text',...
	'String','Analysis',...
	'Position',[025 315 350 020].*WS,...
	'ForegroundColor','k',...
    'FontAngle','Italic',...
	'HorizontalAlignment','Left',...    
	'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 167 360 82].*WS,...
    'Tag', 'EEG')
uicontrol(Fmenu,'Style','Text',...
    'String','Display',...
	'Position',[025 228 350 020].*WS,...
	'ForegroundColor','k',...
    'FontAngle','Italic',...
	'HorizontalAlignment','Left',...
	'Tag', 'EEG')

uicontrol(Fmenu,'Style','Frame',...
	'Position',[020 122 360 40].*WS,...
    'Tag', 'EEG')
uicontrol(Fmenu,'Style','Text',...
	'String','Inference',...
	'Position',[025 142 350 020].*WS,...
	'ForegroundColor','k',...
	'FontAngle','Italic',...
	'HorizontalAlignment','Left',...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','Text',...
	'String','SPM Intracerebral EEG Toolbox',...
	'ToolTipString','modality & defaults set for EEG/MEG',...
	'ForegroundColor',[1 1 1]*.6,'BackgroundColor',[1 1 1]*.8,...
	'FontAngle','Italic','FontWeight','Normal','FontSize',FS(10),...
	'HorizontalAlignment','center',...
	'Position',[020 89 360 020].*WS,...
	'Tag','EEG','Visible','on')

% Utilities frames and text
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','Frame',...
	'BackgroundColor',spm('Colour'),...
	'Position',[010 012 380 75].*WS)

uicontrol(Fmenu,'Style','Text','String',SPMc,...
	'ToolTipString',SPMc,...
	'ForegroundColor',[1 1 1]*.6,...
	'BackgroundColor',[1 1 1]*.8,...
	... 'FontName',PF.times,...
    'FontSize', FS(6), ...
	'HorizontalAlignment','center',...
	'Position',[020 002 360 010].*WS)

%-Objects with Callbacks - main spm_*_ui.m routines
%=======================================================================

%-Data
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','PopUp',...
    'String','Import|  Convert|  Electrodes positions|  Bipolar montage|  Set bad channels|  Detect bad channels|  Events|  Set time origin',...
	'Position',[035 379 100 030].*WS,...
    'FontSize', FS(9), ...
	'ToolTipString', 'Convert data in SPM EEG format',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_spm_eeg_converteeg2mat;','ImaGIN_Electrode;','ImaGIN_BipolarMontage;','ImaGIN_BadChannelSet;','ImaGIN_BadChannel;','ImaGIN_Events;','ImaGIN_TimeZero;'},...
    'Tag', 'EEG');

uicontrol(Fmenu,'Style','PopUp',...
    'String','Preprocess|  Bandpass filter|  Notch filter|  Downsample|  Stimulation artifact: Detect|  Stimulation artifact: Correct',...
	'Position',[150 379 100 030].*WS,... 
    'FontSize', FS(9), ...
	'ToolTipString','Other preprocessing functions',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_BandPassFilter;','ImaGIN_NotchFilter;','ImaGIN_spm_eeg_downsample;','ImaGIN_StimDetect;','ImaGIN_ArtefactCorrectionModel;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','PopUp',...
    'String','Display|  ImaGIN viewer|  SPM viewer|  Brainstorm viewer',...
	'Position',[265 379 100 030].*WS,... 
    'FontSize', FS(9), ...
	'ToolTipString','Other preprocessing functions',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_DispData;','spm_eeg_review;','ImaGIN_DispDataBst;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','Epoch',...
    'Position',[035 349 100 027].*WS,... 
	'ToolTipString','Select a time window according to pre-defined events and overwrite a new (smaller) data set',...
    'UserData','ImaGIN_Crop',...
	'CallBack','ImaGIN_Crop;',...
    'Tag', 'EEG',...
    'Visible','on');

uicontrol(Fmenu,'String','Implantation',...
	'Position',[150 349 100 027].*WS,...
	'ToolTipString','Prepare montage from postop MRI',...
	'UserData','ImaGIN_MRImplantation',...
	'CallBack','ImaGIN_MRImplantation;')

uicontrol(Fmenu,'String','3D interpolation',...
    'Position',[265 349 100 027].*WS,...
	'ToolTipString', '3D interpolation of EEG intracerebral data',...
	'UserData','ImaGIN_spm_eeg_TF_images_3D;',...
	'CallBack','ImaGIN_spm_eeg_TF_images_3D;',...
	'Visible','on',...
    'Tag','EEG')

% Analysis
%-----------------------------------------------------------------------
uicontrol(Fmenu,'Style','PopUp',...
    'String','Time-frequency|  Compute|  Normalise|  Rescale|  Average|FFT',...
	'Position',[035 290 100 030].*WS,...
    'FontSize', FS(9), ...
	'ToolTipString','Time frequency analysis (power and synchrony)',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_spm_eeg_tf;','ImaGIN_NormaliseTF;','ImaGIN_RescaleTF;','ImaGIN_AverageTF;','ImaGIN_FFT;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'Style','PopUp',...
    'String','Causality|  Compute|  Average',...
	'Position',[150 290 100 030].*WS,...
    'FontSize', FS(9), ...
	'ToolTipString','Causality analysis',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'ImaGIN_Causality_compute;','ImaGIN_Causality_average;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','Firing rate',...
	'Position',[035 261 100 027].*WS,...
	'ToolTipString','Firing rate',...
    'CallBack','ImaGIN_FiringRate;',...
	'UserData','ImaGIN_FiringRate',...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','Seizure detect',...
	'Position',[150 261 100 027].*WS,...
	'ToolTipString','Detect seizure and write a time series in a new channel',...
	'UserData','ImaGIN_SeizureDetect',...
	'CallBack','ImaGIN_SeizureDetect;')

uicontrol(Fmenu,'String','Epileptogenicity',...
	'Position',[265 261 100 027].*WS,...
	'ToolTipString','Compute epileptogenicity on SEEG',...
	'UserData','ImaGIN_Epileptogenicity',...
	'CallBack','ImaGIN_Epileptogenicity;')


% Statistical EEG
%-----------------------------------------------------------------------
uicontrol(Fmenu,'String','Data',...
	'Position',[035 205 100 027].*WS,...
	'ToolTipString','Display data',...
	'UserData','ImaGIN_DispData',...
	'CallBack','ImaGIN_DispData;',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Electrodes',...
	'Position',[150 205 100 027].*WS,...
	'ToolTipString','Display positions of electrodes in a SPM .mat/.dat file',...
	'UserData','ImaGIN_DispElectrodes',...
	'CallBack','ImaGIN_DispElectrodes;',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Brain overlay',...
    'Position',[265 205 100 027].*WS,...
	'ToolTipString', '3D EEG intracerebral data overlaid on a brain',...
	'UserData','ImaGIN_Disp3D;',...
	'CallBack','ImaGIN_Disp3D;',...
	'Visible','on',...
    'Tag','EEG')

uicontrol(Fmenu,'String','Time-frequency',...
	'Position',[035 172 100 027].*WS,...
	'ToolTipString','Display time frequency analysis',...
	'UserData','ImaGIN_DispTF',...
	'CallBack','ImaGIN_DispTF',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Causality',...
	'Position',[150 172 100 027].*WS,...
	'ToolTipString','Display causality analysis',...
	'UserData','ImaGIN_DispCausality',...
	'CallBack','ImaGIN_DispCausality',...
    'Tag', 'EEG');

uicontrol(Fmenu,'String','Surface overlay',...
    'Position',[265 172 100 027].*WS,...
	'ToolTipString', 'Display surface results on a cortex surface',...
	'UserData','ImaGIN_DispSurface;',...
	'CallBack','ImaGIN_DispSurface;',...
	'Visible','on',...
    'Tag','EEG')

uicontrol(Fmenu,'String','Results',...
	'Position',[130 127 130 027].*WS,...
    'FontSize',FS(9), ...
	'ToolTipString','Inference and regional responses etc.',...
	'UserData','spm_results_ui',...
	'CallBack','[hReg,xSPM,SPM] = spm_results_ui;',...
    'Tag', 'EEG');

% -Utility buttons (first line)
% -----------------------------------------------------------------------
uicontrol(Fmenu,'Style','PopUp',...
    'String','Display...|images|M/EEG',...
	'Position',[020 054 082 024].*WS,...
	'ToolTipString','orthogonal sections and M/EEG display',...
    'CallBack','spm(''PopUpCB'',gcbo)',...
	'UserData',{'spm_image;','spm_eeg_review;'},...
    'Tag', 'EEG')

uicontrol(Fmenu,'String','Check Reg',...
	'Position',[112 054 083 024].*WS,...
	'ToolTipString','check image registration',...
	'UserData','spm_check_registration',...
	'CallBack','spm_check_registration;')

uicontrol(Fmenu,'Style','PopUp',...
	'String',Modalities,...
	'ToolTipString','change modality PET<->fMRI<->EEG',...
	'Tag','Modality',...
	'Position',[298 054 082 024].*WS,...
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
uicontrol(Fmenu,'String','Quit',...
	'Position',[298 020 082 024].*WS,...
	'ToolTipString','exit SPM',...
	'ForeGroundColor','r',...
	'Interruptible','off',...
	'UserData','QuitSPM',...
	'CallBack','spm(''Quit''), clear all')

end




