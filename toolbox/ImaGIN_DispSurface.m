function ImaGIN_DispSurface()
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
% Authors: Francois Tadel, 2017

%% ===== SELECT DATA =====
% Select surface
SurfaceFile = spm_select(1, '\.gii$', 'Select surface mesh (.gii)');
% Select texture
TextureFile = spm_select(1, '\.gii$', 'Select data texture (.gii)');
% Plot electrodes
tmp = spm_input('Plot electrodes ','+1','Yes|No');
if strcmp(tmp,'Yes')
    FileElec = spm_select(inf, '.*\.(mat|txt)$', 'Select electrode positions file (.txt or .mat/.dat)');
else
    FlagElec = [];
end

%% ===== THRESHOLD =====
% Load texture
gii = gifti(TextureFile);
W = gii.cdata(:,:,:);
% Get minimum/maximum
ma = max(W(:));
mi = min(W(:));
% If we have only positive or negative value: already know what to do
if (ma <= 0)
    threshSign = 'Negative';
elseif (mi >= 0)
    threshSign = 'Positive';
else
    threshSign = spm_input('Sign of the values', '+1', 'Positive|Negative|Both');
end
% Get maximum amplitude
switch (threshSign)
    case 'Negative',  strMax = sprintf('min=%0.3f', mi);
    case 'Positive',  strMax = sprintf('max=%0.3f', ma);
    case 'Both',      strMax = sprintf('max=%0.3f', max(abs(ma), abs(mi)));
end
% Threshold
thresh = str2num(spm_input(['Amplitude threshold (' strMax ') '], '+1', 's'));
if isempty(thresh)
    thresh = 0;
end
% Apply thresold to the data
switch (threshSign)
    case 'Negative'
        W(W >= thresh) = NaN;
    case 'Positive'
        % When displaying delay maps, the background is NaN and 0 means t=0, so we need to include the zero in the threshold
        if (thresh == 0) && (nnz(isnan(W(:))) > 0.2 * numel(W))
            W(W < 0) = NaN;
        % Otherwise: regular thresholding
        else
            W(W <= thresh) = NaN;
        end
    case 'Both'
        W((W >= -thresh) & (W <= thresh)) = NaN;
end


%% ===== DISPLAY SURFACE =====
% Display the mesh and texture
H = spm_mesh_render('Disp', SurfaceFile);
spm_mesh_render('Overlay', H, W);
spm_mesh_render('ColourMap', H, jet);
drawnow;
% Update the light position
camlight(H.light);
% If there are electrodes displayed: make the surface transparent
if ~isempty(FileElec)
    set(H.patch, 'FaceAlpha', 0.7);
end
% Set figure name
set(H.figure, 'Name', TextureFile, 'NumberTitle', 'off');


%% ===== DISPLAY ELECTRODES =====
for i0=1:size(FileElec,1)
    % Get file extension
    [tmp,tmp,fExt] = fileparts(FileElec);
    % Load file
    switch (fExt)
        case '.txt'
            PosElec=load(deblank(FileElec(i0,:)));
        case '.mat'
            D = spm_eeg_load(deblank(FileElec(i0,:)));
            PosElec = struct(D).sensors.eeg.elecpos;
    end
    % Ask for color
    col = spm_input(sprintf('Color for electrode set %d ',i0),'+1','b|g|r|c|m|y|k');
    % Plot all the electrodes
    hold on
    line(PosElec(:,1), PosElec(:,2), PosElec(:,3), ...
         'LineStyle',  'none', ...
         'Marker',     '.', ...
         'MarkerSize', 10, ...
         'Color',      col, ...
         'Parent',     H.axis);
    hold off
end



