function D = ImaGIN_TimeZero(S)
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

try
    Filename = S.Fname;
catch
    Filename = spm_select(1, '\.mat$', 'Select data file');
end
if isempty(Filename)
    D = [];
    return;
end

D = spm_eeg_load(Filename);
Events = events(D);

if iscell(Events)
    Events = Events{:};
end

TimeOnset = timeonset(D);

try
    EventRef = S.EventRef;
catch
    EventRef = spm_input('Name of reference event ', 1, 's');
end
n=0;

for i1=1:length(Events)
    if strcmp(Events(i1).type,EventRef)
        n = n+1;
        Eventref(n) = Events(i1).time;
    end
end

try
    EventRef = min(Eventref);
catch 
    EventRef = 0;
end

try
    Offset = S.Offset;
catch
    Offset = spm_input('Duration [sec] before reference', '+1', 'r',0);
end

D = timeonset(D,-(EventRef-Offset)+TimeOnset);

for i0=1:length(Events)
    Events(i0).time = Events(i0).time-(EventRef-Offset);
end
if isfield(D,'spike')
    for i1 = 1:length(D.spike.timings)
        D.spike.timings{i1} = D.spike.timings{i1}-(EventRef-Offset);
    end
end

% This assigns these events to the first trials (the only one if you have continuous data)
D = events(D, 1, Events);
save(D);

% Add log entry
ImaGIN_save_log(fullfile(D), sprintf('TimeZero: Time offset %fs', -(EventRef-Offset)));



