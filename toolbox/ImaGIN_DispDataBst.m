function ImaGIN_DispDataBst(FileName)
%IMAGIN_DISPDATABST Review .mat/.dat file with Brainstorm
%
% USAGE:  ImaGIN_DispDataBst(FileName)
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
% Authors: Francois Tadel, 2017

% Check Brainstorm installation
if ~exist('brainstorm.m', 'file')
    % Display installation message
    h = msgbox(['Brainstorm is not installed on this computer.' 10 10 ...
            'The download will start automatically in a moment:' 10 ...
            '- Save the brainstorm.zip file on your computer and unzip it,' 10 ...
            '- Add the brainstorm3 folder to your MATLAB path (WITHOUT the subfolders).'], 'Download Brainstorm');
    waitfor(h);
    % Start download
    web('http://neuroimage.usc.edu/bst/getupdate.php?c=UbsM09', '-browser');
    return;
end

% Select file is not specified in input
if (nargin < 1) || isempty(FileName) || exist(FileName, 'file')
    FileName = spm_select(1, '\.mat$', 'Select data file');
    if isempty(FileName)
        return;
    end
end

% Output file name (if modified within Brainstorm)
OutputFile = strrep(FileName, '.mat', '_bst.mat');
% Call Brainstorm to display the file
OutputFile = brainstorm('autopilot', 'ReviewRaw', FileName, 'SPM-DAT', 1, OutputFile);



