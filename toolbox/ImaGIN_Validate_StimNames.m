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

patientCode = S.patientName;

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
    
    xpr4b  = '\w*55.0hz\w*';
    xpr5b  = '\w*55hz\w*';
    xpr6b  = '\w*55 hz\w*';
    
    xpr7  = '\w*alarme\w*';
    xpr8  = '\w*SE1Hz\w*';
    xpr9  = '\w*SE 1Hz\w*';
    xpr10  = 'crise';
    [~,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4)) || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6)) || ~isempty(regexpi(Notes{c},xpr7)) || ...
                ~isempty(regexpi(Notes{c},xpr8)) || ~isempty(regexpi(Notes{c},xpr9)) || ...
                ~isempty(regexpi(Notes{c},xpr4b)) || ~isempty(regexpi(Notes{c},xpr5b)) || ~isempty(regexpi(Notes{c},xpr6b)) ||...
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
    noteName = strrep(noteName,'stim','');
    Notes{KeepEvent(j)} = noteName;
end

rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s'; % find available pulse duration values 

mil_flag = 0;
% MIL patients have stims parameters in stim_parameters file saved on /02-raw
% that could be loaded 
if ~isempty(strfind(patientCode,'MIL'))
    %load('/gin/data/database/02-raw/stim_parameters-ftract-mil.mat','stim_params')
    load('stim_parameters-ftract-mil.mat','stim_params')
    Loc = find(ismember(stim_params.PCode, patientCode), 1);
    if ~isempty(Loc)
        mil_flag = 1;
        Frq = stim_params.Freq{Loc};
        Amp = stim_params.Ampl{Loc};
        Pul = stim_params.Pulse{Loc};
        for n = 1:length(KeepEvent)
            xsub1 = regexp(Notes{KeepEvent(n)},rxp1,'match');
            xsub2 = regexp(Notes{KeepEvent(n)},rxp2,'match');
            xsub3 = regexp(Notes{KeepEvent(n)},rxp3,'match');
            
            if ~isempty(xsub1) % insert stim Amplitude
                Notes{KeepEvent(n)} = char(strrep(Notes{KeepEvent(n)},xsub1,Amp));
            else
                Notes{KeepEvent(n)} = [Notes{KeepEvent(n)} '_' Amp];
            end
            
            if ~isempty(xsub2) % insert stim Frequency
                Notes{KeepEvent(n)} = char(strrep(Notes{KeepEvent(n)},xsub2,Frq));
            else
                Notes{KeepEvent(n)} = [Notes{KeepEvent(n)} '_' Frq];
            end
            
            if ~isempty(xsub3) % insert stim Pulse
                Notes{KeepEvent(n)} = char(strrep(Notes{KeepEvent(n)},xsub3,Pul));
            else
                Notes{KeepEvent(n)} = [Notes{KeepEvent(n)} '_' Pul];
            end
            Notes{KeepEvent(n)} = char(strrep(Notes{KeepEvent(n)},'_3_','_'));
            evt(n).type = Notes{KeepEvent(n)};
        end
        D = events(D,1,evt);
        D2 = clone(D, D.fnamedat, [D.nchannels D.nsamples D.ntrials]);
        D2(:,:,:) = D(:,:,:);
        save(D2);
        fprintf('\n \n MESSAGE: .. %s parameters updated ..::\n',patientCode); 
        set_final_status('OK') 
    end
end
if mil_flag == 0
    pVals = [];
    for k = 1:length(KeepEvent)
        xsub3 = regexp(Notes{KeepEvent(k)},rxp3,'match');
        xsub4 = regexp(Notes{KeepEvent(k)},rxp4,'match');
        if ~isempty(xsub3)&& isempty(xsub4)
            xsub3 = char(xsub3);
            pVals = [pVals;str2double(xsub3(1:end-2))];
        elseif isempty(xsub3)&& ~isempty(xsub4)
            xsub4 = char(xsub4);
            xsub3 = char(xsub4(1:end-1)); % sec --> msec
            if numel(xsub3) <= 3
                buff = num2str(str2double(xsub3)*1000);
                if numel(buff) <= 4
                    xsub3 = buff;
                end
            end
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
            error('No pulse duration found in Notes, no default value neither.');
        else
            if numel(pval) > 1
                pval = pval(1);
            end
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
                Notes{ii} = strrep(Notes{ii},'Stim_Start_',''); Notes{ii} = strrep(Notes{ii},'_0us','');
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
end







