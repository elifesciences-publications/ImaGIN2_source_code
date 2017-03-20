function ImaGIN_deploy()
% IMAGIN_DEPLOY - ImaGIN deployment script.
%
% STEPS:
%    - Update *.m inital comments (with doc/autocomment.txt)
%    - Remove *.asv files
%    - Copy files to GIT folder and open GitGUI

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
% Authors: Francois Tadel, 2017


%% ===== CONFIGURATION =====
% Start timer
tic;
% License file
srcDir = fileparts(fileparts(which(mfilename)));
commentFile = fullfile(srcDir, 'deploy', 'autocomment.txt');


%% ===== GET ALL DIRECTORIES =====
% Get all the subdirectories
srcPath = GetPath(fullfile(srcDir, 'toolbox'));
srcPath = [srcPath, srcDir, ';', fullfile(srcDir, 'deploy')];
% Split string
jPath = java.lang.String(srcPath);
jSplitPath = jPath.split(';');


%% ===== PROCESS DIRECTORIES =====
disp(['DEPLOY> Reading: ', strrep(commentFile, srcDir, '')]);
% Read file
autoComment = ReadAsciiFile(commentFile);
if isempty(autoComment)
    error('Auto-comment file not found.');
end
% Convert to Unix-like string
autoComment = strrep(autoComment, char([13 10]), char(10));
% Initialize line counts
nFiles   = 0;
nCode    = 0;
nComment = 0;
% Updating the M-files
disp('DEPLOY> Updating: Comments in all *.m files...');

for iPath = 1:length(jSplitPath)
    curPath = char(jSplitPath(iPath));
    % Remove ASV files
    delete(fullfile(curPath, '*.asv'));
    
    % === PROCESS M-FILES (EDIT COMMENT, COUNT LINES) ===
    % List all .m files in current directory
    mFiles = dir(fullfile(curPath, '*.m'));
    % Process each m-file
    for iFile = 1:length(mFiles)
        % Build full file name
        fName = fullfile(curPath, mFiles(iFile).name);
        % Replace comment block in file
        ReplaceBlock(fName, '% -==================', '==================-', autoComment);
        % Count files and lines
        [tmpComment, tmpCode] = CountLines(fName, autoComment);
        nFiles   = nFiles + 1;
        nCode    = nCode + tmpComment;
        nComment = nComment + tmpCode;
    end
end
disp('DEPLOY> Statistics:');
disp(['DEPLOY>     Number of files  : ' num2str(nFiles)]);
disp(['DEPLOY>     Lines of code    : ' num2str(nCode)]);
disp(['DEPLOY>     Lines of comment : ' num2str(nComment)]);

% Done
stopTime = toc;
if (stopTime > 60)
    disp(sprintf('DEPLOY> Done in %dmin\n', round(stopTime/60)));
else
    disp(sprintf('DEPLOY> Done in %ds\n', round(stopTime)));
end


% %% ===== COPY TO GIT FOLDER =====
% % Copy all the subfolders
% disp('DEPLOY> Copying to GIT folder...');
% !xcopy C:\Work\Dev\brainstorm3\brainstorm.m C:\Work\Dev\brainstorm_git\brainstorm3\brainstorm.m /y /q
% !xcopy C:\Work\Dev\brainstorm3\defaults C:\Work\Dev\brainstorm_git\brainstorm3\defaults /s /e /y /q
% !xcopy C:\Work\Dev\brainstorm3\deploy   C:\Work\Dev\brainstorm_git\brainstorm3\deploy   /s /e /y /q
% !xcopy C:\Work\Dev\brainstorm3\doc      C:\Work\Dev\brainstorm_git\brainstorm3\doc      /s /e /y /q
% !xcopy C:\Work\Dev\brainstorm3\external C:\Work\Dev\brainstorm_git\brainstorm3\external /s /e /y /q
% !xcopy C:\Work\Dev\brainstorm3\java     C:\Work\Dev\brainstorm_git\brainstorm3\java     /s /e /y /q
% !xcopy C:\Work\Dev\brainstorm3\toolbox  C:\Work\Dev\brainstorm_git\brainstorm3\toolbox  /s /e /y /q
% % Start GIT GUI in the deployment folder
% system('start /b cmd /c ""C:\Program Files\Git\cmd\git-gui.exe" --working-dir "C:\Work\Dev\brainstorm_git\brainstorm3""');



end




%% =================================================================================================
%  ===== HELPER FUNCTIONS ==========================================================================
%  =================================================================================================

%% ===== READ ASCII FILE =====
function fContents = ReadAsciiFile(filename)
    fContents = '';
    % Open ascii file
    fid = fopen(filename, 'r');
    if (fid < 0)
        return;
    end
    % Read file
    fContents = char(fread(fid, Inf, 'char')');
    % Close file
    fclose(fid);
end

%% ===== WRITE ASCII FILE =====
function writeAsciiFile(filename, fContents)
    % Open ascii file
    fid = fopen(filename, 'w');
    if (fid < 0)
        return;
    end
    % Write file
    fwrite(fid, fContents, 'char');
    % Close file
    fclose(fid);
end

%% ===== GET PATH =====
function p = GetPath(d)
    % Generate path based on given root directory
    files = dir(d);
    if isempty(files)
        return
    end
    % Base path: input dir
    p = [d pathsep];
    % Set logical vector for subdirectory entries in d
    isdir = logical(cat(1,files.isdir));
    % Recursively descend through directories
    dirs = files(isdir); % select only directory entries from the current listing
    for i=1:length(dirs)
       dirname = dirs(i).name;
       % Ignore directories starting with '.' and 'defaults' folder
       if (dirname(1) ~= '.') && ~strcmpi(dirname, 'defaults')
           p = [p GetPath(fullfile(d, dirname))]; % recursive calling of this function.
       end
    end
end


%% ===== REPLACE BLOCK IN FILE =====
function ReplaceBlock(fName, strStart, strStop, strNew)
    % Read file
    fContents = ReadAsciiFile(fName);
    % Detect block markers (strStart, strStop)
    % => Start
    iStart = strfind(fContents, strStart);
    if isempty(iStart)
        disp(['*** Block not found in file: "', fName, '"']);
        return;
    end
    iStart = iStart(1);
    % => Stop
    iStop = strfind(fContents(iStart:end), strStop) + length(strStop) + iStart - 1;
    if isempty(iStop)
        disp(['*** Block not found in file: "', strrep(fName,'\','\\'), '"']);
        return;
    end
    iStop = iStop(1);

    % Replace file block with new one
    fContents = [fContents(1:iStart - 1), ...
        strNew, ...
        fContents(iStop:end)];
    % Re-write file
    writeAsciiFile(fName, fContents);
end


%% ===== COUNT LINES =====
function [nComment, nCode] = CountLines(fName, strExclude)
    nComment = 0;
    nCode = 0;
    % Read file
    fContents = ReadAsciiFile(fName);
    % Remove the CR characters
    fContents(fContents == 13) = [];
    % Remove header comment block
    fContents = strrep(fContents, strExclude, '');    
    % Split in lines
    fSplit = str_split(fContents, 10);
    % Loop on lines
    for i = 1:1:length(fSplit)
        fLine = strtrim(fSplit{i});
        if (length(fLine) < 3)
            % Skip
        elseif (fLine(1) == '%')
            nComment = nComment + 1;
        else
            nCode = nCode + 1;
        end
    end
end


