function D = ImaGIN_BadChannel(S)

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
% Authors: Viateur Tuyisenge & Olivier David

% Extract badchannels indices
FileIn = S.dataset; % name of cropped iEEG file .dat/.mat
badDir = S.DirFileOut; % directory where iEEG file .dat/.mat with badchannels indices will be stored
trainDir = S.trainBase;  % directory where training set is stored
try
    D = spm_eeg_load(FileIn); % Load the cropped meeg object
catch
    FileIn = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(FileIn);
end

[~,cutName,~] = fileparts(FileIn);
bPrefix = strcat(badDir,'/',cutName); 

clear S2;
S2.FileName = FileIn;

T = ImaGIN_FeatureSEEG(S2); % function returns a table T of features


Tbase = readtable(strcat(trainDir, '/trainBaseFeatures.csv')); % load training base 
trainingData = Tbase(:,2:9);
[trainedClassifier, ~] = ImaGIN_trainClassifier(trainingData); % train the model 
%Uncomment the following line to use the trained model
%load(strcat(trainDir, '/ImaGIN_trainedClassifier.mat')) %  load the trained classifier
channelClass = trainedClassifier.predictFcn(T(:,2:8)); % predict new dataset
bd = strcmp(channelClass,'Bad');
bIdx = find(bd);
badFile = fopen(strcat(bPrefix,'_bIdx.txt'),'w'); % Save badchannel indices in .txt file
fprintf(badFile,'%d\n',bIdx(:));
fclose(badFile);

Tnew = [T channelClass]; 
Tnew.Properties.VariableNames{'Var9'} = 'Note';
csvfilename = strcat(bPrefix,'.csv'); % Save feature table & badchannels indices
writetable(Tnew,csvfilename,'Delimiter',',');
if ~isempty(bIdx)
    D = badchannels(D,bIdx,1); %Add badchannel index in meeg object
end
Dbad = clone(D, bPrefix, [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
Dbad(:,:,:) = D(:,:,:);
save(Dbad);
%Channel plots and ScreenShots
figDir = strcat(badDir, '/ScreenShot');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

close all;

Size = 8;  %Number of channels per screenshot
tmp = floor(size(D,1)/Size);
for i2 = 1:tmp
    figure(i2);
    set(gcf,'Position',[629 -17 702 1101])
    for i3 = 1:Size
        if intersect(i3+(i2-1)*Size,bIdx) == i3+(i2-1)*Size
            color = 'r'; %Bad channels will be printed in red 
        else
            color = 'k'; %Good channels will be printed in black
        end
        subplot(Size,1,i3)
        plot(time(D),D(i3+(i2-1)*Size,:),color);
        ylabel([num2str(i3+(i2-1)*Size) ' : ' D.chanlabels{i3+(i2-1)*Size}])
        if i3 == 1
            figName = char(strcat(cutName,'_',num2str(i3+(i2-1)*Size),'-', ...
                num2str(i2*Size)));
            title(figName,'interpreter','none');
        end
        axis tight
    end
    zoom on
    fig = figure(i2);
    print(fig,fullfile(figDir,figName),'-dpng'); %ScreenShot
    close;
end

rmd = size(D,1) - tmp*Size;
if rmd ~= 0
    figure(tmp + 1)
    set(gcf,'Position',[629 -17 702 1101])
    for i4 = 1:rmd
        if intersect(i3+(i2-1)*Size+i4,bIdx)== i3+(i2-1)*Size+i4
            color = 'r';
        else
            color = 'k';
        end
        subplot(rmd,1,i4)
        plot(time(D),D(i3+(i2-1)*Size+i4,:),color);
        ylabel([num2str(i3+(i2-1)*Size+i4) ' : ' D.chanlabels{i3+(i2-1)*Size+i4}])
        if i4 == 1
            figName = char(strcat(cutName,'_',num2str(i3+(i2-1)*Size + 1),'-',num2str(size(D,1))));
            title(figName,'interpreter','none');
        end
        axis tight
    end
    zoom on
    fig = figure(i2+1);
    print(fig, fullfile(figDir,figName), '-dpng');
    close
end
disp('Badchannels completed !\n');
