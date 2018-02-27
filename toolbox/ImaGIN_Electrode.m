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
% Copyright (c) 2000-2018 Inserm U1216
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
            % If there are two lists of names
            if iscell(Name{1})
                iChanPos = findChannel(Sensors.label{i1}, Name{1});
                if isempty(iChanPos)
                    iChanPos = findChannel(Sensors.label{i1}, Name{2});
                end
            % If there is only one list of names
            else
                iChanPos = findChannel(Sensors.label{i1}, Name);
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
    
    
    D = sensors(D,'EEG',Sensors);
    
    if ~isempty(S.FileOut)
        D2 = clone(D,S.FileOut,[D.nchannels D.nsamples D.ntrials]); % Create a new .mat/dat file (F-TRACT convention)
        D2(:,:,:) = D(:,:,:);
        save(D2); 
        
        if ~isempty(chFound)
            ImaGIN_save_log(fullfile(D2), 'Positions added for channels:', chFound);
        end
        if ~isempty(chNotFound)
            ImaGIN_save_log(fullfile(D2), 'Positions not found for channels:', chNotFound);
        end
        
    else
        save(D); % Update .mat file        
        % Add entries in log file
        if ~isempty(chFound)
            ImaGIN_save_log(fullfile(D), 'Positions added for channels:', chFound);
        end
        if ~isempty(chNotFound)
            ImaGIN_save_log(fullfile(D), 'Positions not found for channels:', chNotFound);
        end
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


%% Find channel name in a list
function iChanPos = findChannel(Label, List)
    % Look for channel in position file
    iChanPos = find(strcmpi(Label, List));
    % Channel not found: try replacing ' by p
    if isempty(iChanPos)
        iChanPos = find(strcmpi(strrep(Label,'''','p'), strrep(List,'''','p')));
    end
end

