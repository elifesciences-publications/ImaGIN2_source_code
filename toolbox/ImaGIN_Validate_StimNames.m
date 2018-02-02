function ImaGIN_Validate_StimNames(S)
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
% Authors: Viateur Tuyisenge & Olivier David
 
sFile = S.dataset;
pulseDefault = str2double(S.defaultPulseDuration);
fprintf('MESSAGE: Got default pulse duration = %d \n', pulseDefault);
try
    D = spm_eeg_load(sFile); % Load the converted file .mat
catch
    error('File %s NOT loaded \n', sFile);
    %sFile = spm_select(1, '\.mat$', 'Select data file');
    %D=spm_eeg_load(sFile);
end
evt = events(D);
evsize = size(evt,2);        % Number of events
Notes  = cell(1,evsize);     % Events labels
% Extract events properties (label and time in sampling)
for i = 1: evsize
    Notes{i}= char(evt(i).type);
end
% Read notes to keep only those related to a stimulation
KeepEvent=[];
for c=1:evsize % Navigate all available events
    xpr1  = '\w*hz_\w*';
    xpr2  = '\w*stim\w*';
    xpr3  = '\w*mA\w*';
    xpr4  = '\w*50.0hz\w*';
    xpr5  = '\w*50hz\w*';
    xpr6  = '\w*50 hz\w*';
    xpr7  = '\w*alarme\w*';
    xpr8  = '\w*SE1Hz\w*';
    xpr9  = '\w*SE 1Hz\w*';
    xpr10  = 'crise';
    [ds,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4)) || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6)) || ~isempty(regexpi(Notes{c},xpr7)) || ...
                ~isempty(regexpi(Notes{c},xpr8)) || ~isempty(regexpi(Notes{c},xpr9)) || ...
                strcmpi(Notes{c}(1:min([length(Notes{c}) 5])),xpr10)
        elseif ~isempty(regexpi(Notes{c},xpr1))
            KeepEvent=[KeepEvent c];
        elseif ~isempty(regexpi(Notes{c},xpr2))
            KeepEvent=[KeepEvent c];
        elseif ~isempty(regexp(Notes{c},xpr3,'ONCE'))
            KeepEvent=[KeepEvent c];
        elseif ismember(lower(Notes{c}(1)),['a':'z']) && di(1)<=4 && ~strcmp(regexprep(Notes{c},' ',''),'SE1Hz') && ~strcmp(regexprep(Notes{c},' ',''),'SE50Hz')
            KeepEvent=[KeepEvent c];
        end
    end
end
for j=1:length(KeepEvent) % Navigate all stim events
    noteName = strrep(char(Notes{KeepEvent(j)}), ' ','_');
    noteName = regexprep(noteName,'�','u'); %OD
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='';
    noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    
    noteName = strrep(noteName,'.0','');
    noteName = strrep(noteName,'.',''); noteName = strrep(noteName,',','');
    noteName = strrep(noteName,'sec','s');  noteName = strrep(noteName,'AA','A');
    keepN = ''; noteName = strrep(noteName,'stim','');
    try
        fundc = strfind(noteName,'_');
        lNumb = strfind(noteName,noteName(1:fundc(1)-1));
        keepN = noteName(1:fundc(1)-1);
        if(numel(lNumb)) == 2 && ~strcmp(keepN,'A') && ~strcmp(keepN,'H')
            noteName = strrep(noteName,keepN,'CHNAME');
        end
    end
    Notes{KeepEvent(j)} = noteName;
end

rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s'; % find available pulse duration values 
pVals = [];
for k = 1:length(KeepEvent)
    xsub3 = regexp(Notes{KeepEvent(k)},rxp3,'match');
    xsub4 = regexp(Notes{KeepEvent(k)},rxp4,'match');
    if ~isempty(xsub3)&& isempty(xsub4)
        xsub3 = char(xsub3);
        pVals = [pVals;str2double(xsub3(1:end-2))];
    elseif isempty(xsub3)&& ~isempty(xsub4)
        xsub4 = char(xsub4);
        xsub3 = char(xsub4(1:end-1));
        pVals = [pVals;str2double(xsub3)];
    else
        xsub3 = '0';
        pVals = [pVals;str2double(xsub3)];
    end
end
pIdx = find(pVals);

if numel(pIdx) == numel(KeepEvent)
    fprintf('MESSAGE: All Notes include pulse duration ! \n');
    set_final_status('OK')
else
    
    if ~isempty(pIdx)
        pval = unique(pVals(pIdx));
        fprintf('MESSAGE: Pulse duration found, but not in all Notes. Its unique value is %d \n', pval);
    else
        fprintf('MESSAGE: Pulse duration not found in any of the Notes. Using default value = %d \n', pulseDefault);
        pval = pulseDefault;
    end
    
    if isempty(pval)
        error('No pulse duration found in Notes, no default value either.');
    else
        pval = strcat(num2str(pval),'us');
        for c = 1:length(KeepEvent)
            if isempty(strfind(Notes{KeepEvent(c)},pval))
                Notes{KeepEvent(c)}=[Notes{KeepEvent(c)} '_' pval]; % add pulse duration
            end
        end
        nevt = 1;
        for ii =1:evsize
            Notes{ii} = strrep(Notes{ii},'MA','mA');
            Notes{ii} = strrep(Notes{ii},'ma','mA');
            Notes{ii} = strrep(Notes{ii},'Ma','mA');
            Notes{ii} = strrep(Notes{ii},'HZ','Hz');
            Notes{ii} = strrep(Notes{ii},'hz','Hz');
            Notes{ii} = strrep(Notes{ii},'hZ','Hz');
            Notes{ii} = strrep(Notes{ii},'1_HA','1Hz');
            Notes{ii} = strrep(Notes{ii},'_mA','mA');
            Notes{ii} = strrep(Notes{ii},'_Hz','Hz');
            Notes{ii} = strrep(Notes{ii},'Stim_Start_','');
            evt(ii).type = Notes{ii};
            if ~strcmpi(evt(ii).type,'RESET ON') && ~strcmpi(evt(ii).type,'RESET OFF')...
                    && isempty(strfind(evt(ii).type,'Stim_Stop_'))... 
                           && isempty(strfind(evt(ii).type,'OFF'))
                    
                 newEvt(nevt) = evt(ii);
                 nevt = nevt+1;
            end
        end
        D = events(D,1,newEvt);
        D2 = clone(D, D.fnamedat, [D.nchannels D.nsamples D.ntrials]);
        D2(:,:,:) = D(:,:,:);
        save(D2);
        fprintf('\n \n ');
        fprintf('MESSAGE: .. Pulse duration of %s added ..::\n',pval); 
        set_final_status('OK')
    end
end







