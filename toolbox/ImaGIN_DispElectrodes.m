function ImaGIN_DispElectrodes(S)
% Display electrodes positions.
%
% INPUTS:
%    - S: Structure with optional fields
%      |- Fname : Full path to the .mat/.dat file to display
%      |- P     : Full path to the T1 overlay to use in the display

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
% Authors: Olivier David, Francois Tadel


% ===== PARSE INPUTS =====
% Get file to display
if (nargin >= 1) && isfield(S, 'Fname') && ~isempty(S.Fname)
    t = S.Fname;
else
    t = spm_select(1, '\.mat$', 'Select SEEG data file');
    if isempty(t)
        return;
    end
end
% T1 overlay
if (nargin >= 1) && isfield(S, 'Fname') && ~isempty(S.Fname)
    P = S.P;
else
    P = [];
end


% ===== LOAD SEEG POSITIONS =====
% Load file
D = spm_eeg_load(t);
% Get electrodes names and labels
Sensors = sensors(D,'EEG');
% Check for missing information
if isempty(Sensors) || isempty(Sensors.elecpos)
    error('No electrode positions in this file.');
end
% Exclude NaN and (0,0,0) positions 
iGood = find(~any(isnan(Sensors.elecpos),2) & ~all(Sensors.elecpos == 0,2));
if isempty(iGood)
    error('The electrodes positions have not been defined in this file.');
end
% Get positions and labels of the good channels
chLoc   = Sensors.elecpos(iGood,:);
chNames = Sensors.label(iGood);


% ===== GET VOLUME OVERLAY =====
if isempty(P)
    % Template or individual subject?
    tmp = spm_input('T1 overlay ', '+1', 'Template|Subject');
    % Template
    if strcmp(tmp,'Template')
        % Ask for type of atlas to use
        Species = spm_input('Species ','+1','Human|Rat|Mouse');
        % Select corresponding atlas file
        switch Species
            case 'Human'
                P = fullfile(spm('dir'),'canonical','avg152T1.nii');
            case 'Rat'
                P = fullfile(spm_str_manip(spm('dir'),'h'),'ratlas5','template','template_T1.img');
            case 'Mouse' 
                P = fullfile(spm_str_manip(spm('dir'),'h'),'mousatlas5','template','template_T1.img');
        end
    % Individual brain
    else
        P = spm_select(1, 'image', 'Select image for rendering on');
    end
end


% ===== DISPLAY ELECTRODES =====
Nchannels = size(chLoc,1);
% Prepare input structures
sdip.n_seeds = 1;
sdip.n_dip   = Nchannels;
sdip.Mtb     = 1;
sdip.j{1}    = zeros(3 * Nchannels, 1);
sdip.loc{1}  = chLoc;
sdip.Names   = chNames;
% Render electrodes
ImaGIN_DrawElectrodes('Init', sdip, P);


