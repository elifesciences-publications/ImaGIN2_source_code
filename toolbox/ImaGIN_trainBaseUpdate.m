function ImaGIN_trainBaseUpdate(S)

badDir= S.DirIn;

Tbase = readtable('/gin/data/database/02-raw/trainBaseFeatures.csv');
%%
cd(badDir);   % Badchannel directory with csv files (features)
csvTables = dir('*.csv'); % list csv files
for i = 1:length(csvTables)
    csvName = csvTables(i).name;
    csvPath = strcat(badDir,'/',csvName);
    crTable = readtable(csvPath);
    crTable.Properties.VariableNames = {'rankIdx','rankXcorr', 'rankVal', 'ch_dev', 'ch_ampl', 'ch_grad', 'ch_kurt', 'ch_hurs','Note'};
    Tbase = [Tbase;crTable]; 
end

% update trainBase %  /gin/data/database/02-raw/trainBaseFeatures.csv
% trainedModel /gin/data/database/02-raw/ImaGIN_trainedClassifier.mat
trainedfname = '/gin/data/database/02-raw/ImaGIN_trainedClassifier.mat';
csvfname     = '/gin/data/database/02-raw/trainBaseFeatures.csv';
%%
try
    [trainedClassifier, validationAccuracy] = ImaGIN_trainClassifier(Tbase);
    save(trainedfname,'trainedClassifier');
    writetable(Tbase,csvfname,'Delimiter',',');
catch
    writetable(Tbase,csvfname,'Delimiter',',');
end
