function D = ImaGIN_DecoupageAuto(S)

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

sFile = S.dataset;
DirOut= S.DirFileOut;
thisN = S.StimName; % stim event to crop individually 
try
    D = spm_eeg_load(sFile); % Load the converted file .mat
catch
    sFile = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(sFile);
end
%%
evt = events(D);
evsize = size(evt,2);        % Number of events
Notes  = cell(1,evsize);     % Events labels
Time   = zeros(1,evsize);
totTime= max(time(D));
Time0  = min(time(D));
% Extract events properties (label and time in sampling)
for i = 1: evsize
    Notes{i}  = evt(i).type;
    Time(1,i) = evt(i).time;
end
%%
nHz  = 5;  % Accepted stim frequency
minStim = 3; 
bgnTime = zeros(1,evsize); % Beginning of event
endTime = zeros(1,evsize); % End of event

%%
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
    xpr11  = '\w*OFF\w*';
    [ds,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4)) || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6)) || ~isempty(regexpi(Notes{c},xpr7)) || ...
                ~isempty(regexpi(Notes{c},xpr8)) || ~isempty(regexpi(Notes{c},xpr9)) || ...
                ~isempty(regexpi(Notes{c},xpr11)) || strcmpi(Notes{c}(1:min([length(Notes{c}) 5])),xpr10)
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


%% ------------------------------------------------
% Case only specific stim events are to be cropped
if ~isempty(thisN)
    KeepEvent = thisN;
end
%% ------------------------------------------------
for c=1:length(KeepEvent) % Navigate all stim events
    
    % Detect stimulations and stimulation indices & save in text file
    [pth,matFile,~]= fileparts(sFile);
    clear S
    S.Fname = fullfile(pth, matFile);
    S.Channels = [];
    S.StimStart= Time(1,KeepEvent(c))-0.5;
    if c < length(KeepEvent)
        S.StimEnd  = min([Time(1,KeepEvent(c))+180 Time(1,KeepEvent(c+1))-0.5]);
    else
        S.StimEnd  = min([totTime Time(1,KeepEvent(c))+120]);
    end
    S.StimContinuous= 0;
    S.EvtName  = Notes{KeepEvent(c)};
    noteName = strrep(char(Notes{KeepEvent(c)}), ' ','_');
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='';
    noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    noteName = regexprep(noteName,'�','u'); %OD
    noteName = regexprep(noteName,'MA','mA'); %OD
    noteName = regexprep(noteName,'MS','mA'); %OD - errors in Milan notes
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
    %% check if stim electr numbers are concatenated without space or -
    numbr = regexp(noteName,'\d*','Match');
    if numel(numbr) >= 2
        if str2double(numbr(2))~= str2double(numbr(1)) + 1 && str2double(numbr(1))~= str2double(numbr(2)) + 1 
            if numel(numbr{1}) == 2
                elecno = strcat(numbr{1}(1),'-',numbr{1}(2));
                noteName = strrep(noteName,numbr{1},elecno);
            elseif numel(numbr{1}) == 3
                elecno = strcat(numbr{1}(1),'-',numbr{1}(2:3));
                noteName = strrep(noteName,numbr{1},elecno);
            elseif numel(numbr{1}) == 4
                elecno = strcat(numbr{1}(1:2),'-',numbr{1}(3:4));
                noteName = strrep(noteName,numbr{1},elecno);
            end
        end
    end
  %%   
    numb = regexp(noteName,'\d*','Match');
    
    if numel(numb) >= 2
        idx1 = strfind(noteName,numb{1});
        idx2 = strfind(noteName,numb{2});
        cnbre1= numel(numb{1});
        cnbre2= numel(numb{2});
        if str2double(numb(1)) == str2double(numb(2)) - 1
            sfix = noteName(1:(idx1(1)+cnbre1)-1);
            noteName = strcat(upper(sfix),noteName(idx2(1):end));
        elseif str2double(numb(1)) == str2double(numb(2)) + 1
            sfix = noteName(1:idx1(1)-1);
            noteName = strcat(upper(sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
        else
            if str2double(numb(1)) > str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat(upper(sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
            elseif str2double(numb(1)) < str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat(upper(sfix),num2str(numb{1}), num2str(numb{2}),noteName(idx2(1)+cnbre2:end));
            end
        end

        xpr1  = '\w*Hz_\w*'; xpr2 = '\w*us_\w*';
        xpri1 = regexpi(noteName,xpr1); xpri2 = regexpi(noteName,xpr2);
        if isempty(xpri1) && isempty(xpri2)
            endpart  = '_1Hz_Xus';
            noteName = strcat(noteName,endpart);
        end
    end
    noteName = strrep(noteName,'-','');  noteName = strrep(noteName,'_mA','mA');
    
    %% build .mat/.dat name
    
    try
        [ds,di] = regexp(noteName,'\d*','Match');
        xsub0 = noteName(1:(di(1)+numel(ds{1})-1));
        rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
        rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
        xsub1 = regexp(noteName,rxp1,'match'); 
        if isempty(xsub1), xsub1 = '0mA'; 
            xsub1 = cellstr(xsub1);
        end
        xsub2 = regexp(noteName,rxp2,'match'); 
        if isempty(xsub2), xsub2 = '0Hz';
            xsub2 = cellstr(xsub2);
            S.StimFreq = 1;
            StimFreq = 1;
        else
            StimFreq=str2double(xsub2{1}(1:end-2));
            S.StimFreq = StimFreq;
        end
        xsub3 = regexp(noteName,rxp3,'match');
        xsub4 = regexp(noteName,rxp4,'match');
        if ~isempty(xsub3)&& isempty(xsub4)
            xsub3 = cellstr(xsub3);
        elseif isempty(xsub3)&& ~isempty(xsub4)
            xsub4 = char(xsub4);
            xsub3 = cellstr(strcat(xsub4(1:end-1),'us'));
        else
            xsub3 = '0us';
        end
        
        FullN = strcat(upper(xsub0),'_',xsub1,'_',xsub2,'_',xsub3);
        FullN = char(unique(FullN));
        if isempty(FullN)
            FullN = noteName;
        end
    catch
        FullN = noteName;
        S.StimFreq = 1;
    end
    
    noteName = strrep(FullN,'.',',');
    
    %%
    numb2 = regexp(noteName,'\d*','Match');
    if numel(numb2) >= 1
        idxn1 = strfind(noteName,numb2{1});          %OD
        subn1 = strrep(noteName(1:idxn1),'_','');
        noteName = char(strcat(subn1,noteName(idxn1+1:end)));
    end
    if keepN
        noteName = strrep(noteName,'CHNAME',keepN);
    end
    ptrn = ',';
    if strncmp(noteName,ptrn,1)
        noteName = char(noteName(2:end));
    end
    z_sc = strfind(noteName,'_');
    if ~isempty(z_sc)
        zStr = char(regexp(noteName(1:z_sc(1)),'0\d*','Match'));
        if numel(char(zStr)) == 4
            if strcmp(zStr(1),'0') && strcmp(zStr(3),'0')
                nonZ = strrep(zStr,'0','');
                zStr = cellstr(zStr);
                noteName = char(strrep(noteName,zStr,nonZ));
            elseif strcmp(zStr(1),'0') && ~strcmp(zStr(3),'0')
                cpyzStr = cellstr(zStr);
                zStr(1) = '';
                nonZ = cellstr(zStr);
                noteName = char(strrep(noteName,cpyzStr,nonZ));
            end
        end
    end

    %%
    [~,tmpdi]=regexp(noteName,'\d*','Match');
    noteNameNew=noteName;
    noteNameNew(1:tmpdi(1)-1)=upper(noteNameNew(1:tmpdi(1)-1));
    noteNameNew = strrep(noteNameNew,'''','p');     %to avoid ' in the name
    S.FileOut=  fullfile(DirOut, strcat(noteNameNew,'.txt'));
    [stimTime,~,StimulationFreqU] = ImaGIN_StimDetect(S);
    disp(KeepEvent(c)), disp(S.EvtName);
    
    if StimulationFreqU >20
    end
    stimFq = StimFreq;
    
    if numel(stimTime) >= minStim
        [~, stimTimeOut, ~] = fileparts(S.FileOut);
        
        if StimulationFreqU>20  %ugly fix to try to detect 50 Hz stimulation
            stimFq = round(StimulationFreqU);
        else        
            stimStep = stimTime(2:end) - stimTime(1:end-1);
            stimFq   = round(median(stimStep)); % median frequency within event : %OD corrected to median
        end

        % --------------------------------------------------
        % Define a file name for a dissected seeg
        rxp2  = '[-+]?(\d*[.])?\d+Hz';
        FileName = strcat(stimTimeOut,'_1.mat');
        stimHz = regexp(FileName,rxp2,'match');
        strFq = strcat(num2str(stimFq),'Hz');
        if ~strcmp(stimHz,strFq);
            FileName = char(strrep(FileName,stimHz,strFq));
        end
        % Check if the file exists i.e repeated stim
        if exist(fullfile(DirOut, FileName), 'file') ~= 2
            myFileName = FileName;
        else
            numstr = dir(fullfile(DirOut,[FileName(1:end-5),'*','.mat']));
            num    = size(numstr,1) + 1;
            myFileName = sprintf('%s%d.mat',FileName(1:end-5),num);
        end
        namePrefix = myFileName(1:end-4);
        % ----------------------------------------
        % Call ImaGIN_Crop
        bgnTime(KeepEvent(c)) = stimTime(1);
        endTime(KeepEvent(c)) = stimTime(end);
        
        % Select only the stimulations
        if stimFq <= nHz && numel(stimTime) >= minStim && numel(numb) > 1 ...
                && ~strcmp(xsub2(1),'50Hz') 
            clear S
            S.Job   = 'Manual';
            S.Fname = fullfile(pth, matFile);
            S.EventStart= Notes{KeepEvent(c)}; % Note n
            S.numbStart = KeepEvent(c); %%%
            S.numbEnd = KeepEvent(c); %%%
            if c > 1
                if endTime(KeepEvent(c-1)) == 0
                    endTime(KeepEvent(c-1)) = Time(KeepEvent(c-1));
                end
                if Time(KeepEvent(c)) - endTime(KeepEvent(c-1))-1 >= 40
                    setStart = 40;
                else
                    setStart = Time(KeepEvent(c)) - endTime(KeepEvent(c-1))-1;
                    if setStart<1
                        setStart=1;
                    end
                end
            elseif c == 1
                if Time(KeepEvent(c))-Time0 >= 40
                    setStart = 40;
                else
                    setStart = Time(KeepEvent(c)) - Time0;
                end
            end
            S.OffsetStart= setStart;
            S.EventEnd   = Notes{KeepEvent(c)}; % Note n
            S.OffsetEnd  = endTime(KeepEvent(c)) - Time(KeepEvent(c)) + 3;
            S.FileOut     = fullfile(DirOut, namePrefix);
            ImaGIN_Crop(S);
        end
    end
          
    
    if stimFq <= nHz && numel(stimTime) >= minStim && numel(numb) > 1 ...
            && ~strcmp(xsub2(1),'50Hz') ...
            
        % Add events
        clear S2
        S2.Fname   = S.FileOut;
        S2.FileOut = S2.Fname;
        S2.Action  = 'Add';
        S2.Nevent  = 1;
        S2.EventName{1} = 'Stim';
        S2.EventFileName{1} = fullfile(DirOut, strcat(stimTimeOut,'.txt')); % txt file, time in seconds;
        ImaGIN_Events(S2);
        
        % Time origin
        clear S2
        S2.Fname   = S.FileOut;
        S2.FileOut = S2.Fname;
        S2.EventRef= 'Stim';
        S2.Offset  = 0;
        D = ImaGIN_TimeZero(S2);
        
    end
end
f = dir(fullfile(DirOut,'*.txt')); %Delete all text file
f = {f.name};
for k=1:numel(f);
    delete(f{k})
end

disp('Done');
