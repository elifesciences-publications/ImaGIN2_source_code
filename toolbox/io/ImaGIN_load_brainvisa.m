function CortexFile = ImaGIN_load_brainvisa(BvDir, OutputDir)
% IMAGIN_LOAD_BRAINVISA: Concatenante multiple .gii files

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

% Find files
MriFile = {file_find(BvDir, 'nobias_*.nii'), file_find(BvDir, 'nobias_*.nii.gz')};
TessLhFile = file_find(BvDir, '*Lhemi*.gii');
TessRhFile = file_find(BvDir, '*Rhemi*.gii');
TessLhMinf = file_find(BvDir, '*Lhemi*.gii.minf');
TessRhMinf = file_find(BvDir, '*Rhemi*.gii.minf');
% TessLhFile = file_find(BvDir, '*Lwhite*.gii');
% TessRhFile = file_find(BvDir, '*Rwhite*.gii');
% TessLhMinf = file_find(BvDir, '*Lwhite*.gii.minf');
% TessRhMinf = file_find(BvDir, '*Rwhite*.gii.minf');
MriFile(cellfun(@isempty, MriFile)) = [];
if isempty(TessLhFile) || isempty(TessRhFile) || isempty(MriFile) || isempty(TessLhMinf) || isempty(TessRhMinf)
    error('Invalid BrainVISA segmentation folder.');
end
GiiFiles = {TessLhFile, TessRhFile};

% Load first file
outTess = export(gifti(GiiFiles{1}), 'patch');
% Concatenate the other surfaces with the firs one
for i = 2:length(GiiFiles)
    % Load surface
    inTess = export(gifti(GiiFiles{i}), 'patch');
    % Concatenate
    offsetVertices   = size(outTess.vertices,1);
    outTess.faces    = [outTess.faces; inTess.faces + offsetVertices];
    outTess.vertices = [outTess.vertices; inTess.vertices];
end

% % ====================================================
% % OPTION #1: Load the transformation to apply from the .nii (Not working...)
% nii = spm_vol(MriFile{1});
% % Add the transformation to the .gii
% Transf = nii.mat;
% % ====================================================

% ====================================================
% OPTION #2: Read the transformation from the .minf files
% Load minf file
fid = fopen(TessLhMinf,'r');
strMinf = fread(fid, [1 Inf], '*char');
fclose(fid);
% Read transformation from minf
eval([strrep(strMinf, ':', ',') ';']);
iTransf = find(strcmpi(attributes, 'transformations'));
Transf = reshape(attributes{iTransf + 1}(1:16), 4, 4)';
% ====================================================

% Swap faces if the transformation if the transformation introduces a flip of dimensions
if (det(Transf(1:3,1:3)) < 0)
    outTess.faces = outTess.faces(:,[2 1 3]);
end
% Apply transformation to the vertices
if ~isempty(Transf)
    outTess.vertices = bsxfun(@plus, Transf(1:3,1:3) * outTess.vertices', Transf(1:3,4))';
end

% Decimate cortex surface (otherwise computation is too long)
nVertices = 15000;
outTess = reducepatch(outTess, nVertices / length(outTess.vertices));
    
% Output filename
[fPath, fBase, fExt] = fileparts(MriFile{1});
CortexFile = fullfile(OutputDir, [strrep(strrep(strrep(fBase, '.nii', ''), '.gz', ''),'nobias_',''), '.gii']);
% Export new file
gii = gifti(outTess);
save(gii, CortexFile, 'Base64Binary');


