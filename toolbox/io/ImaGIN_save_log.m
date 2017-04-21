function LogFile = ImaGIN_save_log(FileName, Comment, ChanLabels)
% IMAGIN_SAVE_CHANLOG Appends a comment and a list of channel names to a text file.
%
% USAGE:  LogFile = ImaGIN_save_chanlog(FileName, Comment, ChanLabels)
%
% INPUT: 
%    - FileName   : Target file name (LogFile will be the same with _log.txt at the end)
%    - Comment    : Comment to add before the list of channels
%    - ChanLabels : List of channels to list in the file
%
% OUTPUT:
%    - LogFile : Name of the logfile that is generated

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

% Create log filename
[fPath, fBase, fExt] = fileparts(FileName);
LogFile = fullfile(fPath, [fBase, '_log.txt']);
% Open file
fid = fopen(LogFile, 'a+');
if (fid < 0)
    disp(['ImaGIN> Error: Cannot open log file "' LogFile '"']);
    LogFile = [];
    return;
end

% Write time and date
if ~isempty(Comment)
    fprintf(fid, '[%s]\n', char(datetime));
end
% Write comment line
if ~isempty(Comment)
    fprintf(fid, '%s\n', Comment);
end
% Write channel labels
if ~isempty(Comment)
    fprintf(fid, '%s ', ChanLabels{:});
    fprintf(fid, '\n');
end
% Add an empty line to separate the entries
fprintf(fid, '\n');

% Close file
fclose(fid);



