function D = ImaGIN_CropAuto(S)
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
DirOut= S.DirFileOut;
thisN = S.StimName; % stim event to crop individually 
try
    D = spm_eeg_load(sFile); % Load the converted file .mat
catch
    sFile = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(sFile);
end
%% Find existing crop files in DirOut
matExist = dir(fullfile(DirOut,'*.mat')); %Delete all cropped text file
matExist = {matExist.name};

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
    
    xpr4b  = '\w*55.0hz\w*';
    xpr5b  = '\w*55hz\w*';
    xpr6b  = '\w*55 hz\w*';
    
    xpr7  = '\w*alarme\w*';
    xpr8  = '\w*SE1Hz\w*';
    xpr9  = '\w*SE 1Hz\w*';
    xpr10  = 'crise';
    xpr11  = '\w*OFF\w*';
    xpr12  = '\w*t65535\w*';
    [ds,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4)) || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6)) || ~isempty(regexpi(Notes{c},xpr7)) || ...
                ~isempty(regexpi(Notes{c},xpr8)) || ~isempty(regexpi(Notes{c},xpr9)) || ...
                ~isempty(regexpi(Notes{c},xpr12))||~isempty(regexpi(Notes{c},xpr4b)) || ...
                ~isempty(regexpi(Notes{c},xpr5b))|| ~isempty(regexpi(Notes{c},xpr6b))||...
                ~isempty(regexpi(Notes{c},xpr11))|| strcmpi(Notes{c}(1:min([length(Notes{c}) 5])),xpr10)
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
    idxN = strcmp(Notes,thisN);
    idxN= find(idxN);
    KeepEvent = find(KeepEvent==idxN);
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
    noteName = regexprep(noteName,'�sec','us');
    noteName = regexprep(noteName,'�s','us');
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='_';
    noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    noteName = strrep(noteName,'usec','us'); 
    noteName = regexprep(noteName,'MA','mA'); %OD
    noteName = regexprep(noteName,'MS','mA'); %OD - errors in Milan notes
    noteName = regexprep(noteName,'Stim_Start_',''); %YUQ notes
    noteName = regexprep(noteName,'Stim_Stop_','');  %YUQ notes
    noteName = strrep(noteName,'.0','');
    noteName = strrep(noteName,'-','_');  noteName = strrep(noteName,'__','_');    
    noteName = strrep(noteName,'.',''); noteName = strrep(noteName,',','');
    noteName = strrep(noteName,'sec','us');  noteName = strrep(noteName,'_us','us');
    noteName = strrep(noteName,'AA','A');
    keepN = ''; noteName = strrep(noteName,'stim','');  noteName = strrep(noteName,'Stim','');
    try
        fundc = strfind(noteName,'_');
        lNumb = strfind(noteName,noteName(1:fundc(1)-1));
        keepN = noteName(1:fundc(1)-1);
        if(numel(lNumb)) == 2 && ~strcmp(keepN,'A') && ~strcmp(keepN,'H') && ~isempty(keepN)            
            noteName = strrep(noteName,keepN,'CHNAME');            
        end  
    end   
    %% Avoid number starting with 0: 01 02 03,...
    numZ = regexp(noteName,'\d*','Match');
    if numel(numZ) >= 2
        noteName =  char(strrep(noteName,numZ(1), num2str(str2double(numZ(1)))));
        noteName =  char(strrep(noteName,numZ(2), num2str(str2double(numZ(2)))));
    end
    %% check if stim electr numbers are concatenated without space or -
    
    numbr = regexp(noteName,'\d*','Match');
    if numel(numbr) >= 2
        if str2double(numbr(2))~= str2double(numbr(1)) + 1 && str2double(numbr(1))~= str2double(numbr(2)) + 1 
            if numel(numbr{1}) == 2
                 if numel(numbr{2}) ~= 3 && numel(numbr{2}) ~= 2 && numel(numbr{2}) ~= 4
                     if str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                         elecno = strcat(numbr{1}(1),'_',numbr{1}(2));
                         noteName = strrep(noteName,numbr{1},elecno);
                     end
                elseif numel(numbr{2}) == 3                      
                    if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1)) 
                        elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                        noteName = strrep(noteName,numbr{2},elecno);
                    end
                 elseif numel(numbr{2}) == 4
                     if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1))
                         elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3:4));
                         noteName = strrep(noteName,numbr{2},elecno);
                     end
                 elseif numel(numbr{2})== 2
                   if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                      elecno = strcat(numbr{2}(1),'_',numbr{2}(2)); 
                      noteName = strrep(noteName,numbr{2},elecno);
                   end
                end
            elseif numel(numbr{1}) == 3
                if str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3))== 10 
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(2:3))== 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 4
                if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3)) == 10 
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(3)) == 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 5
                 if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4),'_',numbr{1}(5));
                    noteName = strrep(noteName,numbr{1},elecno);
                 end
            elseif numel(numbr{1}) == 1 
                if numel(numbr{2}) == 2
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                        noteName = strrep(noteName,numbr{2},elecno);
                    else
                       elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                       noteName = strrep(noteName,numbr{2},elecno);
                    end
                elseif numel(numbr{2}) == 3
                    elecno = strcat(numbr{2}(1),'_',numbr{2}(2:3));
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1)) 
                        noteName = strrep(noteName,numbr{2},elecno);
                    else
                        elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                        if str2double(numbr(1)) + 1 == str2double(numbr{2}(1:2)) || str2double(numbr(1)) == str2double(numbr{2}(1:2)) + 1
                            noteName = strrep(noteName,numbr{2},elecno);
                        end
                    end
                end
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
        FullN = FullN(1,:);
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
    
    noteName = strrep(noteName,'CHNAME',keepN);  
    
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
    
    %{  
    %% uncomment this section for some datasets
    % stimN1 = str2double(regexp(sfix,'\d*','Match'));
    % xstrN2 = strrep(xsub0,num2str(stimN1),'');
    % stimN2 = str2double(regexp(xstrN2,'\d*','Match'));
    % sfix2  = strrep(sfix,num2str(stimN1), num2str(stimN2));
    % selCh1 = find(strcmpi(elec.label,sfix));
    % selCh2 = find(strcmpi(elec.label,sfix2));
    % S.Channels = [selCh1,selCh2]; 
    % S.FindBadChannels=0;
    %}
    [stimTime,~,StimulationFreqU] = ImaGIN_StimDetect(S);
    disp(KeepEvent(c)), disp(S.EvtName);
    stimFq = StimFreq;
    
    % if numel(stimTime) >= seuilHz         %OD
    if numel(stimTime) >= minStim
        
        [~, stimTimeOut, ~] = fileparts(S.FileOut);
        
%         if StimulationFreqU>20  %OD ugly fix to try to  detect 50 Hz stimulation
%             stimFq = round(StimulationFreqU);
%         else        
            stimStep = stimTime(2:end) - stimTime(1:end-1);
            stimFq   = round(median(stimStep)); % median frequency within event : %OD corrected to median
%         end
              
        % --------------------------------------------------
        % Define a file name for a dissected seeg
        rxp2  = '[-+]?(\d*[.])?\d+Hz';
        FileName = strcat(stimTimeOut,'_1.mat');
        stimHz = regexp(FileName,rxp2,'match');
        strFq = strcat(num2str(stimFq),'Hz');
        if ~strcmp(stimHz,strFq) && ~strcmp(strFq,'0Hz')
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
%%
clc;

f = dir(fullfile(DirOut,'*.mat')); %Delete all cropped text file
f = {f.name};
for k=1:numel(f)
    [~, matFileName, ~] = fileparts(f{k});
    nbrSc = strfind(matFileName,'_');
    txtFileName = matFileName(1:nbrSc(end)-1);
    if exist(fullfile(DirOut,strcat(txtFileName,'.txt')),'file')== 2
        delete(fullfile(DirOut,strcat(txtFileName,'.txt')));
    end
end
nf = dir(fullfile(DirOut,'*.txt'));
nf = {nf.name};

realStims = length(KeepEvent);
realCrps  = numel(f)-numel(matExist);
stimTxt   = numel(nf) + realCrps;
misCrop = 0;
nanStim = 0;

if ~isempty(nf)
    for k2 = 1:numel(nf)
        fID = fopen(fullfile(DirOut, nf{k2}));
        res = {};
        while ~feof(fID)
            res{end+1,1} =fgetl(fID);
        end
        fclose(fID);
        nbrLines = numel(res);
        if nbrLines > minStim
            misCrop = misCrop + 1;
            fprintf('%s not cropped\n',nf{k2});
        else
            nanStim = nanStim +1;
            fprintf('%s is empty \n',nf{k2});
        end
    end
    if misCrop > 0
        fprintf('Number of uncropped Stims:  %d \n',misCrop);
    else
        if nanStim > 0
            fprintf('Crop done completely with %d empty stimulation \n',nanStim);
        else
            disp('Crop done completely');
        end
    end
else
    disp('Crop done completely');
end

fprintf('\n\nNumber of Stims notes found: %d\n',realStims);
fprintf('Number of stimulations detected: %d\n', stimTxt);
fprintf('Number of stimulations cropped: %d\n\n',realCrps);
txtf = dir(fullfile(DirOut,'*.txt')); % delete all .txt file after cropping
txtf = {txtf.name};
for k3 = 1:numel(txtf)
    delete(fullfile(DirOut, txtf{k3}))
end

set_final_status('OK');
disp('Done');

