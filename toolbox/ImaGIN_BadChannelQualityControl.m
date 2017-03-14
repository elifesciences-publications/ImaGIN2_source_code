function D = ImaGIN_BadChannelQualityControl(S)
% FileIn: path linking to meeg object 
FileIn = S.dataset;
try
    D = spm_eeg_load(FileIn); % Load the cropped meeg object
catch
    FileIn = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(FileIn);
end
[badDir,cutName,~] = fileparts(FileIn);

bPrefix = strcat(badDir,'/',cutName); 

try
    bIdx = load(strcat(bPrefix,'_bIdx.txt'));
catch
    bIdx = [];
end

T = readtable(strcat(bPrefix,'.csv'));
tNote = cell(size(T,1),1);
tNote(:) = {'Good'};
tNote(bIdx) = {'Bad'};
T(:,9) = tNote;
csvfilename = strcat(bPrefix,'.csv');
writetable(T,csvfilename,'Delimiter',',');

tmpIdx = D.badchannels; % old badchannels indices
D = badchannels(D,tmpIdx,0); % reset meeg object bad indices

if ~isempty(bIdx)
    D = badchannels(D, bIdx,1); % add new badchannel indeces in meeg object
end
Dbad = clone(D, bPrefix, [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
Dbad(:,:,:) = D(:,:,:);
save(Dbad);

