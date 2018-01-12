function ImaGIN_trainBaseUpdate(badDir, trainBaseFile, trainedFile)
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

% Parse inputs
if (nargin < 3) || isempty(trainedFile)
    trainedFile = '/gin/data/database/02-raw/ImaGIN_trainedClassifier.mat';
end
if (nargin < 2) || isempty(trainBaseFile)
    trainBaseFile = '/gin/data/database/02-raw/trainBaseFeatures.mat';
end

% Read existing train base
if ~isempty(trainBaseFile) && file_exist(trainBaseFile)
    Tbase = load(trainBaseFile);
else
    Tbase.predictors = [];
    Tbase.response = {};
end

% Go to bad channel directory with csv files (features)
cd(badDir);   
% List csv files
csvTables = dir('*.csv');
% Read the csv files and stack them to form the train base
for i = 1:length(csvTables)
    csvName = csvTables(i).name;
    csvPath = strcat(badDir, '/', csvName);
    crTable = readtable(csvPath);
    crTable.Properties.VariableNames = {'rankIdx', 'ch_xcorr', 'ch_var', 'ch_dev', 'ch_ampl', 'ch_grad', 'ch_kurt', 'ch_hurs', 'Note'};
    
    Tbase.predictors = [Tbase.predictors; crTable(:, {'ch_xcorr', 'ch_var', 'ch_dev', 'ch_ampl', 'ch_grad', 'ch_kurt', 'ch_hurs'})];
	Tbase.response   = [Tbase.response; crTable.Note];
end
% Write train base to hard drive
save(trainBaseFile, '-struct', 'Tbase');

% Train classifier
try
    [trainedClassifier, validationAccuracy] = ImaGIN_trainClassifier(Tbase.predictors, Tbase.response);
    save(trainedFile, 'trainedClassifier');
catch
    disp(['ERROR: ' lasterr]);
end


