function D = ImaGIN_BadChannel(S)
% Detect bad channels indices with a trained classifier.
%
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

% Name of input iEEG file .dat/.mat (cropped)
if (nargin >= 1) && isfield(S, 'dataset') && ~isempty(S.dataset)
    FileIn = S.dataset; 
else
    FileIn = spm_select(1, '\.mat$', 'Select EEG mat file');
    if isempty(FileIn)
        return;
    end
end

% Directory where iEEG file .dat/.mat with badchannels indices will be stored
if (nargin >= 1) && isfield(S, 'FileOut') && ~isempty(S.FileOut)
    [badDir, FileOut] = fileparts(S.FileOut);
    isNewFile = 1;
else
    [badDir, FileOut] = fileparts(FileIn);
    isNewFile = 0;
end

% Directory where training set is stored (if not defined, use the toolbox folder)
if (nargin >= 1) && isfield(S, 'trainBase') && ~isempty(S.trainBase)
    trainDir = S.trainBase;
else
    trainDir = fileparts(mfilename('fullpath'));
end


% Detect bad channels and save summary images
% try
    % Load the cropped meeg object
    D = spm_eeg_load(FileIn);
    
    % Compute signal features (apply artifact correction only for stims)
    clear S2;
    S2.FileName = FileIn;
    if ~isempty(D.events) && ~isempty(D.events.type) && ismember('Stim', {D.events.type})
        S2.InterpolationFilter = 1;
    else
        S2.InterpolationFilter = 0;
    end
    T = ImaGIN_FeatureSEEG(S2);
    
    % If the trained classifier is not available: compute it
    trainedFile = fullfile(trainDir, 'ImaGIN_trainedClassifier.mat');
    if ~exist(trainedFile, 'file')
        % Get train base file
        trainBaseFile = fullfile(trainDir, 'ImaGIN_trainBaseFeatures.mat');
        if ~exist(trainBaseFile, 'file')
            error(['File ' trainBaseFile ' was not found.']);
        end
        % Load the train base
        trainBase = load(trainBaseFile);
        % Train the classifier
        trainedClassifier = ImaGIN_trainClassifier(trainBase.predictors, trainBase.response);
        % Save results
        save(trainedFile, '-struct', 'trainedClassifier');
    % Otherwise, load the trained classifier
    else
        trainedClassifier = load(trainedFile);
    end
    
    % Predict new dataset
    channelClass = trainedClassifier.predictFcn(T(:,2:8)); 
    % Get list of detected bad channels
    bIdx = find(strcmp(channelClass, 'Bad'));
       
    %%
    % In case disconnected electrode doesn't have stimulation artefact
    % specific for some FTRACT dataset
    chanLbs = strrep(upper(D.chanlabels) ,'''','p');
    crFname = D.fname;
    undsc   = strfind(crFname,'_');
    if numel(undsc) == 4
        sfix = crFname(1:undsc(1)-1);
        [numb,id] = regexp(sfix,'\d*','Match');
        chl = sfix(1:id(1)-1);
        if numel(numb{1}) == 2
            ch1 = strcat(chl,numb{1}(1));
            ch2 = strcat(chl,numb{1}(2));
        elseif numel(numb{1}) == 3
            ch1 = strcat(chl,numb{1}(1));
            ch2 = strcat(chl,numb{1}(2:3));
        elseif numel(numb{1}) == 4
            ch1 = strcat(chl,numb{1}(1:2));
            ch2 = strcat(chl,numb{1}(3:4));
        end
        chlb = {ch1,ch2};
        chInd  = find(ismember(chanLbs,chlb));
        if ~isempty(chInd)
            if isempty(find(any(bIdx==chInd(1)), 1))
                bIdx(end+1) = [chInd(1)];
            end
            if isempty(find(any(bIdx==chInd(2)), 1))
                bIdx(end+1) = chInd(2);
            end
            bIdx = sort(bIdx);
        end
    end
    %%
    % Save bad channel indices in .txt file
    badFile = fopen(fullfile(badDir, [FileOut, '_bIdx.txt']), 'w'); 
    fprintf(badFile, '%d\n', bIdx(:));
    fclose(badFile);
    
    % Add badchannel index in meeg object
    if ~isempty(bIdx)
        D = badchannels(D,bIdx,1); 
    end
    % Create new file in output (otherwise, simply update the input file)
    if isNewFile
        Dbad = clone(D, fullfile(badDir, FileOut), [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
        Dbad(:,:,:) = D(:,:,:);
        % Save modified list of bad channels
        save(Dbad);
    else
        % Save modified list of bad channels
        save(D);
    end

    % Channel plots and ScreenShots
    figDir = fullfile(badDir, 'ScreenShot');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    
%     close all;
    
    Size = 8;  % Number of channels per screenshot
    n_c  = size(D,1);
    if n_c >= Size
        tmp = floor(n_c/Size);
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
                    figName = char(strcat(FileOut,'_',num2str(i3+(i2-1)*Size),'-', ...
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
                    figName = char(strcat(FileOut,'_',num2str(i3+(i2-1)*Size + 1),'-',num2str(n_c)));
                    title(figName,'interpreter','none');
                end
                axis tight
            end
            zoom on
            fig = figure(i2+1);
            print(fig, fullfile(figDir,figName), '-dpng');
            close
        end
    else
        figure(1)
        set(gcf,'Position',[629 -17 702 1101])
        for i5 = 1:n_c
            if intersect(i5,bIdx)== i5
                color = 'r';
            else
                color = 'k';
            end
            subplot(n_c,1,i5)
            plot(time(D),D(i5,:),color);
            ylabel([num2str(i5) ' : ' D.chanlabels{i5}])
            if i5 == 1
                figName = char(strcat(FileOut,'_',num2str(1),'-',num2str(n_c)));
                title(figName,'interpreter','none');
            end
            axis tight
        end
        zoom on
        fig = figure(1);
        print(fig, fullfile(figDir,figName), '-dpng');
        close
    end
    set_final_status('OK')
    
% catch ME
%     fprintf('%s \n', ME.message)
%     set_final_status('NOK')
% end 
    
