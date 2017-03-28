function ImaGIN_Validate_StimNames(S)
sFile = S.dataset;
clear D; clear D2;
try
    D = spm_eeg_load(sFile); % Load the converted file .mat
catch
    sFile = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(sFile);
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
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='';
    noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    noteName = regexprep(noteName,'�','u'); %OD
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
    
end

rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
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

if ~isempty(pIdx)
    pval = unique(pVals(pIdx));
else
    pval = spm_input('Enter pulse duration',1);
end

clear c;
pval = strcat(num2str(pval),'us');
for c = 1:length(KeepEvent)
    
    [numb,idx] = regexp(Notes{KeepEvent(c)},'\d*','Match');
    
    if numel(numb) > 2
        cnbre1 = numel(numb{1});
        cnbre2 = numel(numb{2});
        if str2double(numb(1)) == str2double(numb(2))-1||str2double(numb(1)) == str2double(numb(2))+1
            if isempty(strfind(Notes{KeepEvent(c)},pval))
                Notes{KeepEvent(c)}=[Notes{KeepEvent(c)} '_' pval]; % add pulse duration
            end
        elseif str2double(numb(1)) < str2double(numb(2))-1
            Notes{KeepEvent(c)}(idx(2):idx(2)+cnbre2-1) = '_';
            Notes{KeepEvent(c)}(idx(2)) =  num2str(str2double(numb(1))+1);
            Notes{KeepEvent(c)}(idx(1):idx(1)+cnbre1-1) = '_';
            Notes{KeepEvent(c)}(idx(1)) =  num2str(str2double(numb(1)));
            if isempty(strfind(Notes{KeepEvent(c)},pval))
                Notes{KeepEvent(c)}=[Notes{KeepEvent(c)} '_' pval ]; % add pulse duration 
            end
        elseif str2double(numb(1)) > str2double(numb(2))+1
            Notes{KeepEvent(c)}(idx(1):idx(1)+cnbre1-1) = '_';
            Notes{KeepEvent(c)}(idx(1)) = num2str(str2double(numb(2))+1);
            Notes{KeepEvent(c)}(idx(2):idx(2)+cnbre2-1) = '_';
            Notes{KeepEvent(c)}(idx(2)) =num2str(str2double(numb(2)));
            if isempty(strfind(Notes{KeepEvent(c)},pval))
                Notes{KeepEvent(c)}=[Notes{KeepEvent(c)} '_' pval]; % add  pulse duration 
            end
        end
    end
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'MA','mA');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'ma','mA');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'Ma','mA');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'HZ','Hz');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'hz','Hz');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'hZ','Hz');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'1_HA','1Hz');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'_mA','mA');
    Notes{KeepEvent(c)} = strrep(Notes{KeepEvent(c)},'_Hz','Hz');
    clear numb; clear cnbre1; clear cnbre2
    evt(KeepEvent(c)).type = Notes{KeepEvent(c)};
end
D = events(D,1,evt);
D2 = clone(D, D.fnamedat, [D.nchannels D.nsamples D.ntrials]);
D2(:,:,:) = D(:,:,:);
save(D2);








