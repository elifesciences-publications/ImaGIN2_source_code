function D = ImaGIN_Electrode(S)
% Set electrode positions.

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

% Get file to edit
try
    t = S.Fname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end
if isempty(t)
    return;
end
P = spm_str_manip(deblank(t(1,:)),'h');

% Read sensor names
try
    Name = S.Name;
catch
    try
        filename = S.filenameName;
        Name = readName(filename);
    catch
        filename = spm_select(1, '\.txt$', 'Select txt file for electrode names', {}, P);
        Name = readName(filename);
    end
end

% Read sensor positions
try
    Position = S.Position;
catch
    try 
        filename = S.filenamePos;
        Position = load(filename);
    catch
        filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
        Position = load(filename);
    end
end

% Set positions
chFound = {};
chNotFound = {};
for i0 = 1:size(t,1)
    T = deblank(t(i0,:));
    D = spm_eeg_load(T);

    Sensors = sensors(D,'EEG');
    Sensors.label = deblank(Sensors.label);
    
    if length(unique(Sensors.label))~=length(Sensors.label)
        warning('repeated label in the TRC file')
    end

    % Loop on all channels available in the file
    for i1 = 1:length(Sensors.chantype)
        if strcmpi(chantype(D,i1),'eeg')
            % Look for channel in position file
            iChanPos = find(strcmpi(Sensors.label{i1}, Name));
            % Channel not found: try replacing ' by p
            if isempty(iChanPos)
                iChanPos = find(strcmpi(strrep(Sensors.label{i1},'''','p'), strrep(Name,'''','p')));
            end
            % If channel was found: get its position
            if ~isempty(iChanPos)
                Sensors.elecpos(i1,:) = Position(iChanPos,:);
                Sensors.chanpos(i1,:) = Position(iChanPos,:);
                chFound{end+1} = Sensors.label{i1};
            else
                disp(['ImaGIN> WARNING: ' Sensors.label{i1} ' not assigned']);
                chNotFound{end+1} = Sensors.label{i1};
            end
        elseif strcmpi(chantype(D,i1),'ecg')
            Sensors.chantype{i1}='ecg';
        end
    end
    
    % Update .mat file
    D = sensors(D,'EEG',Sensors);
    save(D);
    
    % Add entries in log file
    if ~isempty(chFound)
        ImaGIN_save_log(fullfile(D), 'Positions added for channels:', chFound);
    end
    if ~isempty(chNotFound)
        ImaGIN_save_log(fullfile(D), 'Positions not found for channels:', chNotFound);
    end
end

end


%% Read name from file
function Name = readName(filename)
    fid = fopen(filename);
    i1 = 1;
    while 1
        tmp = fgetl(fid);
        if isequal(tmp, -1) || isempty(tmp)
            break;
        end
        Name{i1} = tmp;
        i1 = i1 + 1;
    end
    fclose(fid);
end


