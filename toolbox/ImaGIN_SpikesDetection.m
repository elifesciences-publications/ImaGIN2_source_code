function [spk_Rate] = ImaGIN_SpikesDetection(S)

sFile = S.dataset;
DirOut= S.DirFileOut;
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
StimeWidth = [];
% Extract events properties (label and time in sampling)
for i = 1: evsize
    Notes{i}  = evt(i).type;
    Time(1,i) = evt(i).time;
end

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
    [~,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4)) || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6)) || ~isempty(regexpi(Notes{c},xpr7)) || ...
                ~isempty(regexpi(Notes{c},xpr8)) || ~isempty(regexpi(Notes{c},xpr9)) || ...
                ~isempty(regexpi(Notes{c},xpr12))||~isempty(regexpi(Notes{c},xpr4b)) || ...
                ~isempty(regexpi(Notes{c},xpr5b))|| ~isempty(regexpi(Notes{c},xpr6b))||...
                ~isempty(regexpi(Notes{c},xpr11))|| strcmpi(Notes{c}(1:min([length(Notes{c}) 5])),xpr10)
            KeepEvent=[KeepEvent c];
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
    S.StimFreq = 1; 
    [~,stimIdx,~] = ImaGIN_StimDetect(S); % stims detection
    
    if ~isempty(stimIdx)
        StimeWidth = [StimeWidth;stimIdx(1),stimIdx(end)];
    end
    
    
    disp(KeepEvent(c)), disp(S.EvtName);

end

nonIdx = ~ismember(Notes, Notes(KeepEvent));
nonStimEvents = evt(nonIdx);

idx1 = StimeWidth(1,1) - 3;
SEEG = D(:,1:idx1,1);

for i = 1:size(StimeWidth,1)-1
    wd1 = StimeWidth(i+1,1)-3;
    wd2 = StimeWidth(i,2) + 3;
    SEEG = [SEEG, D(:,wd2:wd1,1)];
end

[~, fname,~] = fileparts(sFile);

basefile  = fullfile(DirOut, ['Baseline_',fname]);
n_t       = size(SEEG,2);
D2        = clone(D, basefile, [D.nchannels n_t 1]);
D2(:,:,1) = SEEG;
D2        = events(D2,1,nonStimEvents);
save(D2);

%%
clearvars -except  basefile;
S.Fname = basefile;
D = ImaGIN_BipolarMontage(S);
% Perform a longitudinal bipolar montage

delete([basefile ,'*'])

%  % TODO: low-pass filter (15 or 30Hz)
%  % TODO: high-pass filter (1.6Hz)
%
%


%%
td = time(D);
idx = find(td);
fs   = D.fsample;
elec   = sensors(D,'eeg');
labels = elec.label;
nonIdx = ~ismember(labels,{'ecg1','ecg2','ECG','ecg'});
labels = labels(nonIdx);
nc     = numel(labels);
seeg   = D(1:nc,(idx));

totRec = (D.nsamples/fs)/60; %time in minutes;

%% Run Delphos_Spikes Detection 
bgnTime = tic;
results = Delphos_detector(seeg,labels, 'SEEG', fs, {'Spk'}, [], [], 40, []);
stpTime = toc(bgnTime);
fprintf('\n Time Delphos spends  %.2f min for  %.2f min recordings ',stpTime/60, totRec);

%% results
labels   = [results.labels]';
spk_Rate = ceil((results.n_Spk)./totRec);
allSPKs  = table(labels, spk_Rate);
%%
pasition = [results.markers.position]';
marks    = {results.markers.label}';
channels = {results.markers.channels}';
%%
SPKmarkers = table(channels,marks, pasition);
[~, bname,~] = fileparts(D.fnamedat);
spkFile = ['spk_' bname];
spkFile.spikeRate = allSPKs;
spkFile.markers   = SPKmarkers;
spkFile.baseline  = ['duration of ' num2str(totRec)  ' minutes'];
spkFile.comment   = 'Spike rate is number of spikes per minute whithin a bipolar channel';

save(fullfile(DirOUt,spkFile), 'spkFile');



