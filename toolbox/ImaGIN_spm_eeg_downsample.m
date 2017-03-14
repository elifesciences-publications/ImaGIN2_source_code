function D = ImaGIN_spm_eeg_downsample(S)
% function used for down-sampling EEG/MEG data
% FORMAT D = spm_eeg_downsample(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file
% Radc_new  - new sampling rate
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_downsample.m 805 2007-04-27 11:26:53Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG downsample setup',0);

try
    DD = S.D;
catch
    DD = spm_select(inf, 'mat', 'Select EEG mat file');
end

for i1=1:size(DD,1)

    D=deblank(DD(i1,:));

    P = spm_str_manip(D, 'H');

    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end

    try
        Radc_new = S.Radc_new;
    catch
        str = 'New sampling rate';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [Radc_new, YPos] = spm_input(str, YPos, 'r');
            if Radc_new < D.fsample, break, end
            str = sprintf('Sampling rate must be less than original (%d)', round(D.fsample));
        end
        S.Radc_new=Radc_new;
    end

    spm('Pointer', 'Watch');drawnow;

    % Prepare for writing data
    
   

    % treat continuous and epoched data differently because of different
    % scaling method

    if size(D, 3) > 1
        % epoched
       % D.scale = zeros(D.Nchannels, 1, D.Nevents);

        spm_progress_bar('Init', D.ntrials, 'Trials downsampled'); drawnow;
        if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials,100));
        else, Ibar = [1:D.ntrials]; end

        nsampl=length(resample(squeeze(D(1, :, 1)), Radc_new, round(D.fsample)));
        try
            Dnew = clone(D, S.FileOut, [D.nchannels nsampl, D.ntrials]);
        catch
            Dnew = clone(D, ['d' fname(D)], [D.nchannels nsampl, D.ntrials]);
        end
        for i2 = 1:D.ntrials
            for i3=1:D.nchannels
            d = squeeze(D(i3, :, i2));
            d2 = resample(d, Radc_new,round(D.fsample));
          %  D.scale(:, 1, i2) = spm_eeg_write(fpd, d2, 2, D.datatype);
            Dnew(i3,:,i2)=d2;
            if ismember(i2, Ibar)
                spm_progress_bar('Set', i2); drawnow;
            end
            end
        end
%         D.events.start = round(D.events.start*Radc_new/D.Radc);
%         D.events.stop = size(d2, 2) - D.events.start - 1;
%         D.Nsamples = size(d2, 2);
    else
        % continuous
      %  D.scale = zeros(D.Nchannels, 1);

        % adjust the timing information
        try
            Events=D.events;
            Events.time = round(Events.time*Radc_new/D.fsample);
        end
        spm_progress_bar('Init', D.nchannels, 'Channels downsampled'); drawnow;
        if D.nchannels > 100, Ibar = floor(linspace(1, D.nchannels, 100));
        else, Ibar = [1:D.nchannels]; end

        nsampl=length(resample(squeeze(D(1,:,:)), Radc_new, round(D.fsample)));
        try
            Dnew = clone(D, S.FileOut, [D.nchannels nsampl, 1]);
        catch
            Dnew = clone(D, ['d' fname(D)], [D.nchannels nsampl, 1]);
        end
        for i = 1:D.nchannels
            d = squeeze(D(i, :, :));
            d2 = resample(d, Radc_new, round(D.fsample));
            Dnew(i,:,1)=d2;
            if ismember(i, Ibar)
                spm_progress_bar('Set', i); drawnow;
            end

        end
       % D.scale = spm_eeg_write(fpd, data, 2, D.datatype);
       
        %ImaGIN
%         D.TimeZero = max([1 round(D.TimeZero*Radc_new/D.Radc)]);
        
    end

    spm_progress_bar('Clear');
    Dnew=fsample(Dnew,Radc_new);
    save(Dnew);

end

spm('Pointer', 'Arrow');
