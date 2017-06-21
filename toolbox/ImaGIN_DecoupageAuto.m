function D = ImaGIN_DecoupageAuto(S)
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
elec   = sensors(D,'EEG');
% Extract events properties (label and time in sampling)
for i = 1: evsize
    Notes{i}  = evt(i).type;
    Time(1,i) = evt(i).time;
end
%%
% Size = 8;  % Number of channels per screenshot
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
    %idxN = strcmp(Notes,thisN);
    %idxN= find(idxN);
    KeepEvent = thisN;%find(KeepEvent==idxN);
end
%% ------------------------------------------------
for c=1:length(KeepEvent) % Navigate all stim events
    
    % Detect stimulations and stimulation indices & save in text file
    [pth,matFile,~]= fileparts(sFile);
    clear S
    S.Filename = fullfile(pth, matFile);
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
    noteName = regexprep(noteName,'MA','mA'); %OD
    noteName = regexprep(noteName,'MS','mA'); %OD - errors in Milan notes
    noteName = strrep(noteName,'.0','');
    noteName = strrep(noteName,'-','_');  noteName = strrep(noteName,'__','_');    
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
        flag = 1;
        idx1 = strfind(noteName,numb{1});
        idx2 = strfind(noteName,numb{2});
        cnbre1= numel(numb{1});
        cnbre2= numel(numb{2});
        %if numel(numb) >= 5
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
            flag = 0;
        %end
        
        if flag == 1
            if str2double(numb(1)) > str2double(numb(2))
                noteName(idx1(1):idx1(1)+cnbre2-1) = numb{2};
                noteName(idx2(1):idx2(1)+cnbre1-1) = numb{1};
            end
            idxs = strfind(noteName(1:idx2(1)),'_');
            if numel(idxs) >= 1
                undsc = strrep(noteName(1:idx2(1)),'_','');
                noteName = strcat(upper(undsc), noteName(idx2(1)+1:end));
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
    
    [~,ix] = regexp(noteName,'\d*','Match');
    if ~isempty(ix)
        alpha  = numel(noteName(1:(ix(1)-1)));
    else
        alpha = numel(noteName);
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
        
        
        % This is already done in ImaGIN_StimDetect
                      %{
        %             % Adjust stims whose note lebel was wrongly placed
        %             % ------------------------------------------------
        %             NoteShift= find(stimStep > nHz);   % Normal stim frequency is 1Hz
        %             if numel(NoteShift) == 1 % but sometimes 2,3,...,10Hz
        %                 mxStim1 = numel(stimTime(1:NoteShift));
        %                 mxStim2 = numel(stimTime(NoteShift+1:end));
        %                 if mxStim1 >= minStim && mxStim2 >= minStim && stimStep(NoteShift)<= 20
        %                     if stimTime(end) < Time(c+1)
        %                         buffer   = stimTime;
        %                         stimTime = buffer;
        %                         stimFq   = round(mean(stimTime(2:NoteShift) - stimTime(1:NoteShift-1)));
        %                         inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
        %                         fprintf(inFile,'%12f\n',stimTime(:));
        %                         fclose(inFile);
        %                     end
        %                 elseif mxStim1 >= minStim && stimStep(NoteShift) > 20
        %                     if stimTime(end) < Time(c+1)
        %                         buffer   = stimTime(1:NoteShift);
        %                         stimTime = buffer;
        %                         stimFq   = round(mean(stimTime(2:end) - stimTime(1:end-1)));
        %                         inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
        %                         fprintf(inFile,'%12f\n',stimTime(:));
        %                         fclose(inFile);
        %                     elseif stimTime(end) >= Time(c+1)
        %                         buffer   = stimTime(1:NoteShift);
        %                         Time(c+1)= stimTime(NoteShift+1)-1;
        %                         stimTime = buffer;
        %                         stimFq   = round(mean(stimTime(2:end) - stimTime(1:end-1)));
        %                         inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
        %                         fprintf(inFile,'%12f\n',stimTime(:));
        %                         fclose(inFile);
        %                     end
        %                 elseif  mxStim2 >= minStim && stimStep(NoteShift)<= 20
        %                     if stimTime(end) < Time(c+1)
        %                         buffer   = stimTime(NoteShift(1)+1:end);
        %                         Time(c)  = stimTime(NoteShift(1)) + 4;
        %                         stimTime = buffer;
        %                         stimFq   = round(mean(stimTime(2:end) - stimTime(1:end-1)));
        %                         inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
        %                         fprintf(inFile,'%12f\n',stimTime(:));
        %                         fclose(inFile);
        %                     end
        %                 end
        %             elseif numel(NoteShift) >= 2
        %                 tmpStim = []; tmpFq = [];
        %                 for ns = 1:numel(NoteShift)
        %                     if stimTime(NoteShift(ns)) < Time(c+1)
        %                         if ns == 1
        %                             mxStim = stimTime(1:NoteShift(1));
        %                             tmpFq1 = mean(mxStim(2:end)-mxStim(1:end-1));
        %                             if numel(mxStim) >= minStim && tmpFq1 <= nHz && stimStep(NoteShift(1)) <= 20
        %                                 tmpFq  = [tmpFq,tmpFq1];
        %                                 tmpStim= [tmpStim,mxStim];
        %                             elseif numel(mxStim) >= 20 && tmpFq1 <= nHz && stimStep(NoteShift(1)) > 20
        %                                 tmpFq  = [tmpFq,tmpFq1];
        %                                 tmpStim= [tmpStim,mxStim];
        %                             end
        %                         elseif ns > 1 && NoteShift(ns)-NoteShift(ns-1) > 3
        %                             mxStim = stimTime(NoteShift(ns-1)+1:NoteShift(ns));
        %                             if mxStim > 1
        %                                 tmpFq1 = mean(mxStim(2:end)-mxStim(1:end-1));
        %                             else
        %                                tmpFq1 = mxStim;
        %                             end
        %                             if numel(mxStim) >= minStim && tmpFq1 <= nHz
        %                                 tmpFq = [tmpFq,tmpFq1];
        %                                 tmpStim = [tmpStim,mxStim];
        %                             end
        %                         end
        %                     elseif stimTime(NoteShift(ns)+1) >= Time(c+1) && tmpFq1 <= nHz
        %                             Time(c+1) = stimTime(NoteShift(ns)+1);
        %                     end
        %                 end
        %                 if numel(tmpStim) >= minStim
        %                     buffer  = tmpStim;
        %                     stimTime= buffer;
        %                     stimFq  = round(mean(tmpFq));
        %                     Time(c) = stimTime(1) - 4;
        %                     inFile  = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
        %                     fprintf(inFile,'%12f\n',stimTime(:));
        %                     fclose(inFile);
        %                 end
        %              end
        
        %}
        
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
                % && ~strcmpi(Notes{c},'PI') ...
                % && isempty(strfind(Notes{c},'PI')) %&&  alpha <= 4 ...
                % && isempty(strfind(Notes{c},'SE'))
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
            %                 %% Rename a cropped meeg object
            %                 DC = spm_eeg_load([S.Prefix matFile]);
            %                 Dnew = clone(DC, S.Prefix, [DC.nchannels DC.nsamples DC.ntrials]);
            %                 Dnew(:,:,:) = DC(:,:,:); save(Dnew);
            %                 delete(strcat(S.Prefix,matFile,'.mat'));
            %                 delete(strcat(S.Prefix,matFile,'.dat'));
        end
    end
          
         %{
    %     else
    %         % Stim detection of the last event
    %         [pth, matFile, ~] = fileparts(sFile);
    %         S.Filename = fullfile(pth, matFile);
    %         S.Channels = [];
    %         S.StimStart= Time(1,evsize) - 4;
    %         S.StimEnd  = totTime;
    %         S.StimFreq = 1;
    %         S.StimContinuous= 0;
    %         S.EvtName  = Notes{evsize};
    %         %%
    %         noteName = strrep(char(Notes{evsize}), ' ','_');
    %
    %         noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='';
    %         noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    %         noteName = strrep(noteName,'.0','');   noteName = strrep(noteName,'''','p');
    %         noteName = strrep(noteName,'.',''); noteName = strrep(noteName,',','');
    %         noteName = strrep(noteName,'sec','s'); noteName = strrep(noteName,'AA','A');
    %         %%
    %         keepN = '';
    %         try
    %             fundc = strfind(noteName,'_');
    %             lNumb = strfind(noteName,noteName(1:fundc(1)-1));
    %             keepN = noteName(1:fundc(1)-1);
    %             if(numel(lNumb)) == 2
    %                 noteName = strrep(noteName,keepN,'CHNAME');
    %             end
    %         end
    %         %%
    %         numb = regexp(noteName,'\d*','Match');
    %
    %         if numel(numb) >= 2
    %             flag = 1;
    %             idx1 = strfind(noteName,numb(1));
    %             idx2 = strfind(noteName,numb(2));
    %             cnbre1= numel(numb{1});
    %             cnbre2= numel(numb{2});
    %             if numel(numb) >= 5
    %                 if str2double(numb(1)) == str2double(numb(2)) - 1
    %                     sfix = noteName(1:(idx1(1)+cnbre1));
    %                     noteName = strcat(upper(sfix),noteName(idx2(1):end));
    %                 elseif str2double(numb(1)) == str2double(numb(2)) + 1
    %                     sfix = noteName(1:idx1(1)-1);
    %                     noteName = strcat(upper(sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
    %                 end
    %                 flag = 0;
    %             end
    %
    %             if flag == 1
    %                 if str2double(numb(1)) > str2double(numb(2))
    %                     noteName(idx1(1)) = numb{2};
    %                     noteName(idx2(1)) = numb{1};
    %                 end
    %                 idxs = strfind(noteName(1:idx2(1)),'_');
    %                 if numel(idxs) >= 1
    %                     undsc = strrep(noteName(1:idx2(1)),'_','');
    %                     noteName = strcat(upper(undsc), noteName(idx2(1)+1:end));
    %                 end
    %             end
    %             xpr1  = '\w*Hz_\w*'; xpr2 = '\w*us_\w*';
    %             xpri1 = regexpi(noteName,xpr1); xpri2 = regexpi(noteName,xpr2);
    %             if isempty(xpri1) && isempty(xpri2)
    %                 endpart  = '_1Hz_Xus';
    %                 noteName = strcat(noteName,endpart);
    %             end
    %         end
    %        noteName = strrep(noteName,'-','');
    %        noteName =  strrep(noteName,'_mA','mA');
    %        noteName = strrep(noteName,'stim','');
    %        %% build .mat/.dat name
    %        try
    %            [ds,di] = regexp(noteName,'\d*','Match');
    %            xsub0 = noteName(1:(di(1)+numel(ds{1})-1));
    %            rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
    %            rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
    %            xsub1 = regexp(noteName,rxp1,'match'); if isempty(xsub1), xsub1 = '0mA'; xsub1 = cellstr(xsub1);end
    %            xsub2 = regexp(noteName,rxp2,'match'); if isempty(xsub2), xsub2 = '0Hz'; xsub2 = cellstr(xsub2);end
    %            xsub3 = regexp(noteName,rxp3,'match');
    %            xsub4 = regexp(noteName,rxp4,'match');
    %            if ~isempty(xsub3)&& isempty(xsub4)
    %                xsub3 = cellstr(xsub3);
    %            elseif isempty(xsub3)&& ~isempty(xsub4)
    %                xsub4 = char(xsub4);
    %                xsub3 = cellstr(strcat(xsub4(1:end-1),'us'));
    %            else
    %                xsub3 = '0us';
    %            end
    %            FullN = strcat(upper(xsub0),'_',xsub1,'_',xsub2,'_',xsub3);
    %            FullN = char(unique(FullN));
    %            if isempty(FullN)
    %                FullN = noteName;
    %            end
    %        catch
    %            FullN = noteName;
    %        end
    %        try
    %            noteName = strrep(FullN,'.',',');
    %        catch
    %            noteName = FullN;
    %        end
    %
    %        numb2 = regexp(noteName,'\d*','Match');
    %        if numel(numb2) >= 1
    %            idxn1 = strfind(noteName,numb2(1));
    %            subn1 = strrep(noteName(1:idxn1),'_','');
    %            noteName = char(strcat(subn1,noteName(idxn1+1:end)));
    %        end
    %        if keepN
    %            noteName = strrep(noteName,'CHNAME',keepN);
    %        end
    %        ptrn = ',';
    %        if strncmp(noteName,ptrn,1)
    %            noteName = char(noteName(2:end));
    %        end
    %        %%
    %        z_sc = strfind(noteName,'_');
    %        if ~isempty(z_sc)
    %            zStr = char(regexp(noteName(1:z_sc(1)),'0\d*','Match'));
    %            if numel(char(zStr)) == 4
    %                if strcmp(zStr(1),'0') && strcmp(zStr(3),'0')
    %                    nonZ = strrep(zStr,'0','');
    %                    zStr = cellstr(zStr);
    %                    noteName = char(strrep(noteName,zStr,nonZ));
    %                elseif strcmp(zStr(1),'0') && ~strcmp(zStr(3),'0')
    %                   cpyzStr = cellstr(zStr);
    %                   zStr(1) = '';
    %                   nonZ = cellstr(zStr);
    %                   noteName = char(strrep(noteName,cpyzStr,nonZ));
    %                end
    %            end
    %        end
    %        [~,ix] = regexp(noteName,'\d*','Match');
    %
    %        if ~isempty(ix)
    %            alpha  = numel(noteName(1:(ix(1)-1)));
    %        else
    %            alpha = numel(noteName);
    %        end
    %         %%
    %         S.FileOut= fullfile(DirOut, strcat(noteName,'.txt'));
    %
    %         [stimTime, ~] = ImaGIN_StimDetect(S);
    %
    %         disp(evsize);
    %
    %         if  numel(stimTime) >= seuilHz
    %              % Define a file name for a dissected seeg
    %             [pnme, stimTimeOut, ~] = fileparts(S.FileOut);
    %             stimStep = stimTime(2:end) - stimTime(1:end-1);
    %             stimFq   = round(mean(stimStep));
    %             NoteShift= find(stimStep > nHz);   % Normal stim frequency is 1Hz
    %             if numel(NoteShift) == 1 % but sometimes 2,3,...,10Hz
    %                 buffer   = stimTime(1:NoteShift(1));
    %                 stimTime = buffer;
    %                 stimFq   = round(mean(stimTime(2:end) - stimTime(1:end-1)));
    %                 inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
    %                 fprintf(inFile,'%12f\n',stimTime(:));
    %                 fclose(inFile);
    %             elseif numel(NoteShift) == 2
    %                 for ns = 1:numel(NoteShift)
    %                     mxStim = numel(stimTime(1:NoteShift(ns)));
    %                     if mxStim >= minStim
    %                         if ns == 2
    %                             stIndx = NoteShift(ns-1);
    %                             endIndx= NoteShift(ns);
    %                         elseif ns == 1
    %                             endIndx= NoteShift(2);
    %                             stIndx = 0;
    %                         end
    %                         %%break;
    %                     else
    %                         endIndx= numel(stimTime);
    %                         stIndx = NoteShift(2);
    %                     end
    %                 end
    %                 buffer = stimTime(stIndx+1:endIndx);
    %                 Time(evsize)= stimTime(stIndx+1) - 2;
    %                 stimTime = buffer;
    %                 stimFq   = round(mean(stimTime(2:end) - stimTime(1:end-1)));
    %                 inFile   = fopen(fullfile(pnme,[stimTimeOut,'.txt']),'w');
    %                 fprintf(inFile,'%12f\n',stimTime(:));
    %                 fclose(inFile);
    %              end
    %             % ----------------
    %
    %             %FileName = strcat(stimTimeOut,'_',num2str(stimFq),'Hz_1.mat');
    %             rxp2  = '[-+]?(\d*[.])?\d+Hz';
    %             FileName = strcat(stimTimeOut,'_1.mat');
    %             stimHz = regexp(FileName,rxp2,'match');
    %             strFq = strcat(num2str(stimFq),'Hz');
    %             if ~strcmp(stimHz,strFq);
    %                  FileName = char(strrep(FileName,stimHz,strFq));
    %             end
    %
    %             % Check if the file exists i.e repeated stim
    %             if exist(fullfile(DirOut,FileName), 'file') ~= 2
    %                myFileName = FileName;
    %             else
    %                 numstr = dir(fullfile(DirOut,[FileName(1:end-5),'*','.mat']));
    %                 num    = size(numstr,1) + 1;
    %                 myFileName = sprintf('%s%d.mat',FileName(1:end-5),num);
    %             end
    %             namePrefix  = myFileName(1:end-4);
    %
    %             % ----------------------------------------
    %             % Call ImaGIN_Crop for the last event
    %             bgnTime(evsize) = stimTime(1);
    %             endTime(evsize) = stimTime(end);
    %
    %             if stimFq <= nHz && numel(stimTime) > minStim && numel(numb) > 1 ...
    %                     && ~strcmpi(Notes{c},'PI') && ~strcmp(xsub2,'50Hz')...
    %                     && isempty(strfind(Notes{c},'PI')) && alpha <= 4 ...
    %                     && isempty(strfind(Notes{c},'SE'))
    %                 clear S
    %                 S.Job   = 'Manual';
    %                 S.Fname = fullfile(pth, matFile);
    %                 S.EventStart = Notes{evsize};  % Note n
    %                 S.numb = c; %%%%
    %                 if evsize > 1
    %                     if Time(evsize) - endTime(evsize-1) >= 40
    %                         setStart = 40;
    %                     else
    %                         setStart = Time(evsize) - endTime(evsize-1) - 2;
    %                     end
    %                 elseif evsize == 1
    %                     if Time(evsize) -Time0 >= 40
    %                         setStart = 40;
    %                     else
    %                         setStart = Time(evsize) - Time0 - 2;
    %                     end
    %                 end
    %                 S.OffsetStart= setStart;
    %                 S.EventEnd   = Notes{evsize};  % End Note
    %                 S.OffsetEnd  = endTime(evsize) - Time(evsize) + 4;
    %                 S.Prefix     = fullfile(DirOut, namePrefix);
    %                 ImaGIN_Crop(S);
    %                 %% Rename a cropped meeg object
    %                 DC = spm_eeg_load([S.Prefix matFile]);
    %                 Dnew = clone(DC, S.Prefix, [DC.nchannels DC.nsamples DC.ntrials]);
    %                 Dnew(:,:,:) = DC(:,:,:); save(Dnew);
    %                 delete(strcat(S.Prefix,matFile,'.mat'));
    %                 delete(strcat(S.Prefix,matFile,'.dat'));
    %
    %             end
    %         end
    %     end
    
    %     if numel(stimTime) > minStim && numel(numb) > 1 && ~strcmpi(Notes{c},'PI') ...
    %             && ~strcmp(xsub2,'50Hz') && isempty(strfind(Notes{c},'PI')) ...
    %             &&  alpha <= 4 && isempty(strfind(Notes{c},'SE'))
    %}
    
    if stimFq <= nHz && numel(stimTime) >= minStim && numel(numb) > 1 ...
            && ~strcmp(xsub2(1),'50Hz') ...
            
        % Add events
        clear S2
        S2.Filename= S.FileOut;
        S2.FileOut = S2.Filename;
        S2.Action  = 'Add';
        S2.Nevent  = 1;
        S2.EventName{1} = 'Stim';
        S2.EventFileName{1} = fullfile(DirOut, strcat(stimTimeOut,'.txt')); % txt file, time in seconds;
        ImaGIN_Events(S2);
        
        % Time origin
        clear S2
        S2.Filename= S.FileOut;
        S2.FileOut = S2.Filename;
        S2.EventRef= 'Stim';
        S2.Offset  = 0;
        D = ImaGIN_TimeZero(S2);
        %{        
%         % File name correction if wrong note (try to detect stimulation
%         % channels)
%         clear S3
%         S3.D = S2.FileOut;
%         S3.ZeroCrossing = 0;
%         S3.Ratio = [];
%         S3.TimeRange = [[-2.5 -0.5]; [min(time(D)) min(time(D))+2]];
%         S3.FrequencyResolution = 2;
%         S3.FrequencyRange = [];
%         D2   = ImaGIN_FFT(S3);
%         elec = sensors(D,'EEG');
%         [figBase,cutName,~] = fileparts(S3.D);
%         %read-out
%         FR=zeros(1,nchannels(D2));
%         for i1=1:nchannels(D2)-1
%             Sp1 = abs(D2.fft.Data{1}(i1,:));
%             Sp2 = abs(D2.fft.Data{2}(i1,:));
%             FR(i1) = Sp1(25)/Sp2(25); % frequences ratio
%         end
%         %Channel of stimulation
%         sfix  = strfind(cutName,'_');
%         numb = regexp(cutName(1:sfix-1),'\d*','Match');
%         chanletter=upper(cutName(1:length(sfix)-length(numb{1})));
%         if length(numb{1})==2
%             Chan1=strcat(chanletter,numb{1}(1));
%             Chan2=strcat(chanletter,numb{1}(2));
%         elseif length(numb{1})==3
%             Chan1=strcat(chanletter,numb{1}(1));
%             Chan2=strcat(chanletter,numb{1}(2:3));
%         elseif length(numb{1})
%             Chan1=strcat(chanletter,numb{1}(1:2));
%             Chan2=strcat(chanletter,numb{1}(3:4));
%         end
%         for i1=1:length(elec.label)
%             if strcmp(Chan1,upper(elec.label{i1}))
%                 FR1=FR(i1);
%                 i1
%             end
%             if strcmp(Chan2,upper(elec.label{i1}))
%                 FR2=FR(i1);
%                 i1
%             end
%         end
%         %Need to find a readout that works to detect stimulation, this one
%         %does not work
%         if ??
%             sfix  = strfind(cutName,'_');
%             suffix= cutName(sfix(1):end);
%             CorrectName = strcat(figBase,'/',chan,chan2(2:end),suffix);
%             % Rename .mat & .dat files;
%             DC = spm_eeg_load(S3.D);
%             Dnew = clone(DC,CorrectName, [DC.nchannels DC.nsamples DC.ntrials]);
%             Dnew(:,:,:) = DC(:,:,:); save(Dnew);
%             delete(strcat(S3.D,'.mat'));
%             delete(strcat(S3.D,'.dat'));
%             disp('File name changed');
%             S2.FileOut=CorrectName;
%         end
         %}        
        
        %% Compute channel features to extract badchannels
        %{  
%         clear S2
%         
%         [cropDir,cutName,~] = fileparts(S.FileOut);
%         shIdx  = strfind(cropDir,'/');
%         badDir = strcat(cropDir(1:shIdx(end)),'BadChannel');
%         if ~exist(badDir, 'dir')
%             mkdir(badDir);
%         end
%         
%         bPrefix = strcat(badDir,'/',cutName); 
%         S2.FileName = S.FileOut;
%         T = ImaGIN_FeatureSEEG(S2); % function returns a table T of features
%         csvfilename = strcat(bPrefix,'.csv'); % Save feature table in BadChannel directory
%         writetable(T,csvfilename,'Delimiter',',');
%                 
%         load ImaGIN_trainedClassifier.mat % path??? load the trained classifier
%         yfit = trainedClassifier.predictFcn(T(:,2:8)); 
%         bd = strcmp(yfit,'Bad');
%         bIdx = find(bd);
%        
%         badFile = fopen(strcat(bPrefix,'_bIdx.txt'),'w');
%         fprintf(badFile,'%d\n',bIdx(:));
%         fclose(badFile); 
        
%         D = badchannels(D,bIdx,1); % add badchannel index in meeg object
%         Dbad = clone(D, bPrefix, [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
%         Dbad(:,:,:) = D(:,:,:);
%         save(Dbad);
%         % Channel plots and ScreenShots 
%         elec = sensors(D,'EEG');
%         figDir = strcat(badDir, '/ScreenShot');
%         if ~exist(figDir, 'dir')
%             mkdir(figDir);
%         end
%         
%         close all;
%         tmp = floor(size(D,1)/Size);
%        
%         for i2 = 1:tmp
%             figure(i2);
%             set(gcf,'Position',[629 -17 702 1101])
%             for i3 = 1:Size
%                 if (i3+(i2-1)*Size) <= size(yfit,1)
%                     if strcmp(yfit(i3+(i2-1)*Size),'Bad')
%                         color = 'r';
%                     else
%                         color = 'k';
%                     end
%                 else
%                     color = 'k';
%                 end
%                 subplot(Size,1,i3)
%                 plot(time(D),D(i3+(i2-1)*Size,:),color);
%                 ylabel([num2str(i3+(i2-1)*Size) ' : ' elec.label{i3+(i2-1)*Size}])
%                 if i3 == 1
%                     figName = char(strcat(cutName,'_',num2str(i3+(i2-1)*Size),'-', ...
%                         num2str(i2*Size)));
%                     title(figName,'interpreter','none');
%                 end
%                 axis tight
%             end
%             zoom on
%             fig = figure(i2);
%             print(fig,fullfile(figDir,figName),'-dpng'); %ScreenShot
%             close;
%         end
%         
%         rmd = size(D,1) - tmp*Size;
%         if rmd ~= 0
%             figure(tmp + 1)
%             set(gcf,'Position',[629 -17 702 1101])
%             for i4 = 1:rmd
%                 if (i3+(i2-1)*Size+i4) <= size(yfit,1)
%                     if strcmp(yfit(i3+(i2-1)*Size+i4),'Bad')
%                         color = 'r';
%                     else
%                         color = 'k';
%                     end
%                 else
%                     color = 'k';
%                 end
%                 subplot(rmd,1,i4)
%                 plot(time(D),D(i3+(i2-1)*Size+i4,:),color);
%                 ylabel([num2str(i3+(i2-1)*Size+i4) ' : ' elec.label{i3+(i2-1)*Size+i4}])
%                 if i4 == 1
%                     figName = char(strcat(cutName,'_',num2str(i3+(i2-1)*Size + 1),'-',num2str(size(D,1))));
%                     title(figName,'interpreter','none');
%                 end
%                 axis tight
%             end
%             zoom on
%             fig = figure(i2+1);
%             print(fig, fullfile(figDir,figName), '-dpng');
%             close
%         end
%         % -----------------------------------------------------------------
       %} 
    end
end
%%
f = dir(fullfile(DirOut,'*.mat')); %Delete all cropped text file
f = {f.name};
for k=1:numel(f);
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

fprintf('Number of existing stimulations:  %d \n',realStims);

if ~isempty(nf)
    misCrop = 0;
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
            fprintf(' %s   not cropped\n',nf{k2});
        end
    end
    if misCrop > 0
        fprintf('Number of uncropped Stims:  %d \n',misCrop);
    else
        disp('Crop done completely');
    end
else
    disp('Crop done completely');
end

txtf = dir(fullfile(DirOut,'*.txt')); % delete all .txt file after cropping
txtf = {txtf.name};
for k3 = 1:numel(txtf)
    delete(fullfile(DirOut, txtf{k3}))
end

