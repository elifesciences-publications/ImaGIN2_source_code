function ImaGIN_Epileptogenicity(S)
% Compute epileptogenicity using time-windowed fft
%
% USAGE:   D = ImaGIN_Epileptogenicity(S)

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
% Authors: Olivier David

warning off
NameEpileptogenicity='EI';

try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select EEG mat file');
end
try
    BB = S.B;
catch
    BB = spm_select(inf, '\.mat$', 'Select Baseline EEG mat file');
end

try
    latency=S.Latency;
catch
    latency = spm_input('Peri-onset time [s]', 1, 'r', 0:4:20, inf);
end
latency=transpose(latency(:));

try
    Atlas = S.Atlas;
catch
    str   = 'Select atlas';
    Atlas = spm_input(str, '+1','Human|Rat|Mouse|PPN');
end

try
    FreqBand=S.FreqBand;
catch
    FreqBand = spm_input('Frequency Band [Hz]', '+1', 'r', [60 100], 2);
end

try
    CorticalMesh = S.CorticalMesh;
catch
    str   = 'Use cortical mesh ';
    str=spm_input(str, '+1','Yes|No');
    if strcmp(str,'Yes')
        CorticalMesh = 1;
    else
        CorticalMesh = 0;
    end
end

if CorticalMesh
    try
        sMRI = S.sMRI;
    catch
        sMRI = spm_select(Inf, 'image', 'Select normalised MRI');
    end
end

try
    Horizon=S.HorizonT;
catch
    Horizon = spm_input('Mesoscopic time scale [s]', '+1', 'r', 4, 1);
end

try
    TimeResolution=S.TimeResolution;
catch
    TimeResolution = spm_input('Time resolution [s]', '+1', 'r', 0.2, 1);
end

try
    ThDelay=S.ThDelay;
catch
    ThDelay = spm_input('Propagation threshold', '+1', 'r', 0.05, 1);
end

try
    AR=S.AR;
catch
    tmp = spm_input('AR correction', '+1', 'Yes|No');
    if strcmp(tmp,'Yes')
        AR=1;
    else
        AR=0;
    end
end

try
    FileName=S.FileName;
catch
    FileName = spm_input('File name', '+1', 's');
end

% %Bad channels
% try
%     BadChannel=S.BadChannel;
% catch
%     BadChannel = spm_input('Bad channels', '+1', 'i','0');
% end
% if sum(BadChannel)==0
%     BadChannel=[];
% end

%Time window
if length(Horizon)==1
    TimeWindow=[0:TimeResolution:Horizon+1+max(latency(:))];
%     Start=min(TimeWindow)-1/0.05;
%     End=max(TimeWindow)+1/0.05;
end
% TimeWindowStart=TimeWindow(1);
Freq50=[48:52 98:102 148:152 198:202 248:252 298:302 348:352 398:402 448:452 498:502];

% %Transform BadChannels into cells
% if ~iscell(BadChannel)
%     tmp=BadChannel;
%     BadChannel={};
%     for i0=1:size(DD,1)
%         BadChannel{i0}=tmp;
%     end
% end
    

%find common Channels and define as bad the missing ones over files
N = zeros(1,size(DD,1));
Labels = cell(1,size(DD,1));
BadChannel = cell(1,size(DD,1));
for i0=1:size(DD,1)
    D=spm_eeg_load(deblank(DD(i0,:)));
    Labels{i0}=chanlabels(D);
    N(i0)=length(Labels{i0});
    BadChannel{i0}=badchannels(D);
end
L=zeros(size(DD,1),max(N));
for i0=1:size(DD,1)
    tmp=setdiff(1:size(DD,1),i0);
    for i1=1:length(Labels{i0})
        for i2=tmp
            for i3=1:length(Labels{i2})
                if strcmp(Labels{i0}(i1),Labels{i2}(i3))
                    L(i0,i1)=L(i0,i1)+1;
                end
            end
        end
    end
end
M=max(L(:));
for i0=1:size(DD,1)
    BadChannel{i0}=unique([BadChannel{i0} find(L(i0,:)<M)]);
    BadChannel{i0}=BadChannel{i0}(find(BadChannel{i0}<=N(i0)));
end

                    
        




for i00=1:size(latency,2)
    
    Latency=mean(latency(:,i00));
    
    for i0=1:size(DD,1)
        
        if length(Horizon)==1
            TimeWindow=[0:TimeResolution:Horizon+1+max(latency(:))];
        elseif length(Horizon)>1
            TimeWindow=[0:TimeResolution:Horizon(i0)+1+max(latency(:))];
%             Start=min(TimeWindow)-1/0.05;
%             End=max(TimeWindow)+1/0.05;
        end
        
        D=spm_eeg_load(deblank(DD(i0,:)));
        Dinit=D;
        P=spm_str_manip(deblank(DD(i0,:)),'h');
        cd(P)
        time8=time(D);
        
%         for i2=1:D.ntrials
%             try
%                 Events=D.events{i2};
%             catch
%                 Events=D.events;
%             end
%             for i1=1:size(Events,2)
%                 if strcmp(Events(i1).type,'Start')
%                     TimeWindowStart=Events(i1).time-30;
%                 end
%                 if strcmp(Events(i1).type,'End')
%                     TimeWindowEnd=DEvents(i1).time;
%                 end
%             end
%         end
                
        %Baseline
        if ~isempty(BB)
            B=spm_eeg_load(deblank(BB(i0,:)));
            Binit=B;
            timebaseline=time(B);
            TimeWindowBaseline=[timebaseline(1):(TimeWindow(2)-TimeWindow(1)):(timebaseline(end)-1)];
        end
        
        %Downsample data in time
        Coarse=1;
        while D.fsample/Coarse>2*max(FreqBand)
            Coarse=Coarse+1;
        end
        Coarse=Coarse-1;
        Data=D(:,1:Coarse:D.nsamples);
        time8=time8(1:Coarse:end);
        if ~isempty(BB)
            databaseline=B(:,1:Coarse:end);
            timebaseline=timebaseline(1:Coarse:end);
        end
        
%         %Crop the data
%         Start=unique(find(abs(time8-Start)==min(abs(time8-Start))));
%         End=unique(find(abs(time8-End)==min(abs(time8-End))));
%         Data=Data(:,Start:End);
%         Time=time8(Start:End);
        
        %compute power using multitaper
        
        clear SS
        SS.D=deblank(DD(i0,:));
        SS.Pre=['Epi_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' FileName];
        SS.Method='Multitaper';
        SS.TimeResolution=TimeResolution;
        SS.frequencies=setdiff(min(FreqBand):max(FreqBand),Freq50);
        SS.frequencies=min(FreqBand):max(FreqBand);
        SS.TimeWindow=[min(TimeWindow) max(TimeWindow)];
        SS.TimeWindowWidth=1;
        SS.channels=1:D.nchannels;
        SS.FactMod=10;
        SS.NSegments=1;
        SS.Taper='hanning';
        try
            DPower=spm_eeg_load(fullfile(D.path,['m1_' SS.Pre '_' D.fname]));
            DPowerNorm=spm_eeg_load(fullfile(D.path,['nm1_' SS.Pre '_' D.fname]));
            if ~isempty(BB)
                DPowerBaseline=spm_eeg_load(fullfile(B.path,['m1_' SS.Pre '_' B.fname]));
            end
        catch
            ImaGIN_spm_eeg_tf(SS);
            DPower=spm_eeg_load(fullfile(D.path,['m1_' SS.Pre '_' D.fname]));
            if isempty(BB)
                SS2.TimeWindow=DPower.tf.time;
                SS2.D=fullfile(D.path,['m1_' SS.Pre '_' D.fname]);
                SS2.B=[TimeWindow(1) TimeWindow(max(find(SS2.TimeWindow<SS2.TimeWindow(1)-SS2.TimeWindow(2))))];
                ImaGIN_NormaliseTF(SS2);
            else
                SSB=SS;
                SSB.D=deblank(BB(i0,:));
                SSB.TimeWindow=TimeWindowBaseline;
                ImaGIN_spm_eeg_tf(SSB);
                DPowerBaseline=spm_eeg_load(fullfile(B.path,['m1_' SSB.Pre '_' B.fname]));
                SS2.D=fullfile(D.path,['m1_' SS.Pre '_' D.fname]);
                SS2.B=fullfile(B.path,['m1_' SSB.Pre '_' B.fname]);
                ImaGIN_NormaliseTF(SS2);
            end
            DPowerNorm=spm_eeg_load(fullfile(D.path,['nm1_' SS.Pre '_' D.fname]));
        end
        
%         clear SS
%         SS.D=deblank(DD(i0,:));
%         SS.Synchro='No';
%         SS.Pre=['Epi_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' FileName];
%         SS.Method='Morlet wavelet';
%         SS.frequencies=setdiff(min(FreqBand):max(FreqBand),Freq50);
%         SS.FactMod=0;
%         SS.Mfactor=40;
%         SS.Width=0;
%         SS.TimeWindow=TimeWindow;
%         SS.TimeWindowWidth=TimeWindow(2)-TimeWindow(1);
%         SS.Coarse=0;
%         SS.channels=1:D.nchannels;
%         SS.TimeResolution=TimeResolutionTF;
%         try
%             DPower=spm_eeg_load(fullfile(D.path,['w1_' SS.Pre '_' D.fname]));
%             DPowerNorm=spm_eeg_load(fullfile(D.path,['nw1_' SS.Pre '_' D.fname]));
%             if ~isempty(BB)
%                 DPowerBaseline=spm_eeg_load(fullfile(B.path,['w1_' SS.Pre '_' B.fname]));
%             end
%         catch
%             ImaGIN_spm_eeg_tf(SS);
%             DPower=spm_eeg_load(fullfile(D.path,['w1_' SS.Pre '_' D.fname]));
%             if isempty(BB)
%                 TimeWindow=DPower.tf.time;
%                 SS2.D=fullfile(D.path,['w1_' SS.Pre '_' D.fname]);
%                 SS2.B=[TimeWindow(1) TimeWindow(max(find(TimeWindow<TimeWindow(1)-TimeWindow(2))))];
%                 ImaGIN_NormaliseTF(SS2);
%             else
%                 SSB=SS;
%                 SSB.D=deblank(BB(i0,:));
%                 SSB.TimeWindow=TimeWindowBaseline;
%                 ImaGIN_spm_eeg_tf(SSB);
%                 DPowerBaseline=spm_eeg_load(fullfile(B.path,['w1_' SSB.Pre '_' B.fname]));
%                 SS2.D=fullfile(D.path,['w1_' SS.Pre '_' D.fname]);
%                 SS2.B=fullfile(B.path,['w1_' SSB.Pre '_' B.fname]);
%                 ImaGIN_NormaliseTF(SS2);
%             end
%             DPowerNorm=spm_eeg_load(fullfile(D.path,['nw1_' SS.Pre '_' D.fname]));
%         end
        
        
        Power=DPower(:,:,:);
        PowerNorm=DPowerNorm(:,:,:);
        PowerBaseline=DPowerBaseline(:,:,:);
        TimeWindow=DPower.tf.time;
        
        %Find frequency band
        IndexFreq1=min(find(SS.frequencies>=min(FreqBand))):max(find(SS.frequencies<=max(FreqBand)));
                
        %Compute power within frequencies of interest
        Epileptogenicity=squeeze(mean(Power(:,IndexFreq1,:),2));
        EpileptogenicityBaseline=squeeze(mean(PowerBaseline(:,IndexFreq1,:),2));
%         if ~isempty(BadChannel)
%             Epileptogenicity(BadChannel(find(BadChannel<=size(Epileptogenicity,1))),:)=NaN;
%             EpileptogenicityBaseline(BadChannel(find(BadChannel<=size(EpileptogenicityBaseline,1))),:)=NaN;
%         end
        if ~isempty(BadChannel{i0})
            Epileptogenicity(BadChannel{i0},:)=NaN;
            EpileptogenicityBaseline(BadChannel{i0},:)=NaN;
        end
        Epileptogenicity=log(Epileptogenicity);
        EpileptogenicityBaseline=log(EpileptogenicityBaseline);
        %Add offset to have only positive values as for fMRI
        %(otherwise problem with globals calculation)
        tmp=min([Epileptogenicity(:);EpileptogenicityBaseline(:)]);
        Epileptogenicity=Epileptogenicity-tmp;
        EpileptogenicityBaseline=EpileptogenicityBaseline-tmp;
                
        %Save Log Power
%         D1=D;
%         D1=path(D,P);
        D1=clone(D,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.mat'],[D.nchannels size(Epileptogenicity,2) 1]);
        D1(:,:,:)=Epileptogenicity;
        D1=fsample(D1,1/(DPower.tf.time(2)-DPower.tf.time(1)));
        D1=timeonset(D1,min(DPower.tf.time));
        save(D1);
%         D1=B;
%         D1=path(B,P);
        D1=clone(B,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.mat'],[B.nchannels size(EpileptogenicityBaseline,2) 1]);
        D1(:,:,:)=EpileptogenicityBaseline;
        D1=fsample(D1,1/(DPowerBaseline.tf.time(2)-DPowerBaseline.tf.time(1)));
        D1=timeonset(D1,min(DPowerBaseline.tf.time));
        save(D1);
        
        %Write 3D images of log power for statistics
        clear SS
        SS.n=3;
        try
            SS.TimeWindow=latency(i0,i00)+[0:TimeResolution:Horizon];
        catch
            SS.TimeWindow=Latency+[0:TimeResolution:Horizon];
        end
        SS.TimeWindowWidth=0;
        SS.interpolate_bad=0;
        SS.SizeSphere=5;
        SS.SizeHorizon=10;
%         SS.SizeHorizon=15;
        SS.CorticalMesh=CorticalMesh;
        if CorticalMesh
            SS.sMRI=sMRI;
        end
        SS.Atlas=Atlas;
        SS.SaveMNI=0;
        if isdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
            cd(P)
            rmdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
        end
        SS.Fname=fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]);
        ImaGIN_spm_eeg_convertmat2ana_3D(SS)
        clear SS
        SS.n=3;
%         SS.TimeWindow=min(DPowerBaseline.tf.time):(TimeResolution*(max(DPowerBaseline.tf.time)-min(DPowerBaseline.tf.time))/(2*Horizon)):max(DPowerBaseline.tf.time);
        SS.TimeWindow=min(DPowerBaseline.tf.time):TimeResolution:max(DPowerBaseline.tf.time);
        SS.TimeWindowWidth=0;
        SS.interpolate_bad=0;
        SS.SizeSphere=5;
        SS.SizeHorizon=10;
%         SS.SizeHorizon=15;
        SS.CorticalMesh=CorticalMesh;
        if CorticalMesh
            SS.sMRI=sMRI;
        end
        SS.Atlas=Atlas;
        SS.SaveMNI=0;
        if isdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
            cd(P)
            rmdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
        end
        SS.Fname=fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]);
        ImaGIN_spm_eeg_convertmat2ana_3D(SS)
        
        %smooth images to get Gaussian fields
        [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
        for i1=1:size(files,1)
            tmp=deblank(files(i1,:));
            Q = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],tmp);
%             M=spm_vol(Q);
%             I=spm_read_vols(M);
%             spm_smooth(M,M.fname,[5 5 5]);
%             M2=spm_vol(Q);
%             I2=spm_read_vols(M2);
%             I2(find(isnan(I)))=NaN;
%             I2=(max(I(:))./max(I2(:)))*I2;
%             spm_write_vol(M2,I2);
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = {[Q ',1']};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 1;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            spm('defaults', 'EEG');
            spm_jobman('run', matlabbatch);
            movefile(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],['s' tmp]),fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],tmp))
        end
        [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
        for i1=1:size(files,1)
            tmp=deblank(files(i1,:));
            Q = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],tmp);
%             M=spm_vol(Q);
%             I=spm_read_vols(M);
%             spm_smooth(M,M.fname,[5 5 5]);
%             M2=spm_vol(Q);
%             I2=spm_read_vols(M2);
%             I2(find(isnan(I)))=NaN;
%             I2=(max(I(:))./max(I2(:)))*I2;
%             spm_write_vol(M2,I2);
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = {[Q ',1']};
            matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 1;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
            spm('defaults', 'EEG');
            spm_jobman('run', matlabbatch);
            movefile(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],['s' tmp]),fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],tmp))
        end

%         %whiten data
%         %baseline
%         [xX] = whiten_data(TimeResolution*(max(DPowerBaseline.tf.time)-min(DPowerBaseline.tf.time))/(2*Horizon),Inf,ones(1,size(files,1)),fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(round(mean(Horizon))) '_' num2str(round(Latency))]),files);
        
        %SPMs
        if isdir(fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
            cd(P)
            rmdir(fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
        end
        mkdir(fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))])};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units='scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT=TimeResolution;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t=16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0=1;
        [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
        ntmp=size(files,1);
        for i1=1:size(files,1)
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i1,1} = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],[files(i1,:) ',1']);
        end
        [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
        for i1=1:size(files,1)
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans{ntmp+i1,1} = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],[files(i1,:) ',1']);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond=struct('name',{},'onset',{},'duration',{},'tmod',{},'pmod',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi{1}='';
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress.name='Seizure';
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress.val=[ones(ntmp,1);zeros(size(files,1),1)];
        matlabbatch{1}.spm.stats.fmri_spec.sess.multireg{1}='';
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf=Inf;
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs=[0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt=1;
        matlabbatch{1}.spm.stats.fmri_spec.global='None';
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        if AR
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        else
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';
        end
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat')};
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat = {fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat')};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '+';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        spm_get_defaults('mask.thresh', 0) %no implicit masking of SPM-Ts
        spm_jobman('run',matlabbatch)
        
        %Write T values for each electrode
        load(fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat'));
        V=spm_vol(fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'spmT_0001.nii'));
        VV=spm_read_vols(V);
        tmp=spm('Defaults','EEG');
        bb=tmp.normalise.write.bb;
        [x,y,z] = meshgrid(bb(1,1):SS.n:bb(2,1),...
            bb(1,2):SS.n:bb(2,2),...
            bb(1,3):SS.n:bb(2,3));
        VV2=permute(VV,[2 1 3]);
        tmp=sensors(D,'EEG');
        try
            PosElec=tmp.elecpos';
        catch
            PosElec=tmp.pnt';
        end
        IndElec=zeros(size(PosElec,2),1);
        EIGamma=zeros(size(PosElec,2),1);
        for i1=1:size(PosElec,2)
            dist=(x(:)-PosElec(1,i1)).^2+(y(:)-PosElec(2,i1)).^2+(z(:)-PosElec(3,i1)).^2;
%             [tmp1,tmp2]=min(dist);
%             IndElec(i1)=tmp2(1);
%             EIGamma(i1)=VV2(IndElec(i1));
            IndElec=find(dist<100);      %interpolate up to 10 mm (SizeHorizon when creating images)
            tmp1=VV2(IndElec);
            EIGamma(i1)=mean(tmp1(tmp1~=0));
        end
%         D1=D;
%         D1=path(D,P);
        D1=clone(D,[NameEpileptogenicity '_' spm_str_manip(Dinit.fname,'s') '_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.mat'],[size(PosElec,2) 1 1]);
        D1(:,:,:)=EIGamma;
        D1=timeonset(D1,0);
        save(D1);
        
        %Write text file
        tmp=fname(D);
        E=spm_eeg_load(deblank(DD(1,:)));
        Montage=sensors(E,'EEG');
        Montage.Nchannels=length(Montage.label);
        SaveFile=fullfile(P,[NameEpileptogenicity '_' tmp(1:end-4) '_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.txt']);
        fid = fopen(SaveFile,'wt');
        fprintf(fid,[' Electrode /  ' NameEpileptogenicity '     \n']);
        fprintf(fid,'\n');
        for i1=1:Montage.Nchannels
            try
                if ~isnan(EIGamma(i1))
                    fprintf(fid,'%s %10.2f \n', cell2mat(Montage.label{i1}), EIGamma(i1));
                else
                    fprintf(fid,'%s NaN \n', cell2mat(Montage.label{i1}));
                end
            catch
                if ~isnan(EIGamma(i1))
                    fprintf(fid,'%s %10.2f \n', cell2mat(Montage.label(i1)), EIGamma(i1));
                else
                    fprintf(fid,'%s NaN \n', cell2mat(Montage.label(i1)));
                end
            end
        end
        fclose(fid);
        
    end
    
    
    if size(DD,1)>1
        
        %Group
        %SPM of GI
        if isdir(fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
            cd(P)
            rmdir(fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
        end
        mkdir(fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]))
        clear matlabbatch
        matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))])};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units='scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT=TimeResolution;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t=16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0=1;
        for i0=1:size(DD,1)
            D=spm_eeg_load(deblank(DD(i0,:)));
            Dinit=D;
            P=spm_str_manip(deblank(DD(i0,:)),'h');
            [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
            ntmp=size(files,1);
            for i1=1:size(files,1)
                matlabbatch{1}.spm.stats.fmri_spec.sess(i0).scans{i1,1} = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],[files(i1,:) ',1']);
            end
            [files,dirs] = spm_select('List',fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'.nii');
            for i1=1:size(files,1)
                matlabbatch{1}.spm.stats.fmri_spec.sess(i0).scans{ntmp+i1,1} = fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],[files(i1,:) ',1']);
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).cond=struct('name',{},'onset',{},'duration',{},'tmod',{},'pmod',{});
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).multi{1}='';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).regress.name='Seizure';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).regress.val=[ones(ntmp,1);zeros(size(files,1),1)];
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).multireg{1}='';
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).hpf=Inf;
            matlabbatch{1}.spm.stats.fmri_spec.sess(i0).fact=struct('name',{},'levels',{});
        end
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs=[0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt=1;
        matlabbatch{1}.spm.stats.fmri_spec.global='None';
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        if AR
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        else
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';
        end
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat')};
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.stats.con.spmmat = {fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat')};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '+';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = ones(1,size(DD,1));
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;
        spm_get_defaults('mask.thresh', 0) %no implicit masking of SPM-Ts
        spm_jobman('run',matlabbatch)
        
        %Write T values for each electrode
        load(fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'SPM.mat'));
        V=spm_vol(fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'spmT_0001.nii'));
        VV=spm_read_vols(V);
        tmp=spm('Defaults','EEG');
        bb=tmp.normalise.write.bb;
        [x,y,z] = meshgrid(bb(1,1):SS.n:bb(2,1),...
            bb(1,2):SS.n:bb(2,2),...
            bb(1,3):SS.n:bb(2,3));
        VV2=permute(VV,[2 1 3]);
        tmp=sensors(D,'EEG');
        try
            PosElec=tmp.elecpos';
        catch
            PosElec=tmp.pnt';
        end
        IndElec=zeros(size(PosElec,2),1);
        EIGamma=zeros(size(PosElec,2),1);
        for i1=1:size(PosElec,2)
            dist=(x(:)-PosElec(1,i1)).^2+(y(:)-PosElec(2,i1)).^2+(z(:)-PosElec(3,i1)).^2;
%             [tmp1,tmp2]=min(dist);
%             IndElec(i1)=tmp2(1);
%             EIGamma(i1)=VV2(IndElec(i1));
            IndElec=find(dist<100);      %interpolate up to 10 mm (SizeHorizon when creating images)
            tmp1=VV2(IndElec);
            EIGamma(i1)=mean(tmp1(tmp1~=0));
        end
%         D1=D;
%         D1=path(D,P);
        D1=clone(D,[NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.mat'],[size(PosElec,2) 1 1]);
        D1(:,:,:)=EIGamma;
        %   D=dtype(D,'float32');
        D1=timeonset(D1,0);
        save(D1);
        
        %Write text file
        E=spm_eeg_load(deblank(DD(1,:)));
        Montage=sensors(E,'EEG');
        Montage.Nchannels=length(Montage.label);
        SaveFile=fullfile(P,[NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency))) '.txt']);
        fid = fopen(SaveFile,'wt');
        fprintf(fid,[' Electrode /  ' NameEpileptogenicity  ' \n']);
        fprintf(fid,'\n');
        for i1=1:Montage.Nchannels
            try
                if ~isnan(EIGamma(i1))
                    fprintf(fid,'%s %10.2f \n', cell2mat(Montage.label{i1}), EIGamma(i1));
                else
                    fprintf(fid,'%s NaN \n', cell2mat(Montage.label{i1}));
                end
            catch
                if ~isnan(EIGamma(i1))
                    fprintf(fid,'%s %10.2f \n', cell2mat(Montage.label(i1)), EIGamma(i1));
                else
                    fprintf(fid,'%s NaN \n', cell2mat(Montage.label(i1)));
                end
            end
        end
        fclose(fid);
        
    end
    
    %delete temporary files
    for i0=1:size(DD,1)
        D=spm_eeg_load(deblank(DD(i0,:)));
        Dinit=D;
        P=spm_str_manip(deblank(DD(i0,:)),'h');
        rmdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
        rmdir(fullfile(P,[FileName spm_str_manip(Dinit.fname,'s') '_' NameEpileptogenicity 'Baseline_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))]),'s')
    end
    
end

% Compute map of propagation delay (only if more than one latency)
if length(latency) > 1

    for i0=1:size(DD,1)
        Delay=NaN*zeros(53,63,46);
        Delay=NaN*zeros(V.dim(1),V.dim(2),V.dim(3));
        D=spm_eeg_load(deblank(DD(i0,:)));
        Dinit=D;
        P=spm_str_manip(deblank(DD(i0,:)),'h');
        for i2=1:size(latency,2)
            Latency=mean(latency(:,i2));
            SPMFile=fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))], 'SPM.mat');
            load(SPMFile)
            % Load spmT map
            P1 = spm_vol(fullfile(P,['SPM_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s') '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(mean(Latency)))],'spmT_0001.nii'));
            Vb1 = spm_read_vols(P1);
            % Calcule threshold for statistical map
            df = [SPM.xCon(1).eidf SPM.xX.erdf];
            S    = SPM.xVol.S;                  %-search Volume {voxels}
            R    = SPM.xVol.R;                  %-search Volume {resels}
    %         u = spm_uc_FDR(0.001,df,'T',1,P1);
            u = spm_uc(ThDelay,df,'T',R,1,S);
            % Activated voxels
            Q1=find(Vb1>=u);
            Q2=find(isnan(Delay));
            Q3=intersect(Q2,Q1);
            Delay(Q3)=mean(latency(:,i2));
        end
        P0=P1;
        P0.fname=fullfile(P,['Delay_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s')  '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']);
        P0 = spm_write_vol(P0,Delay);
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[P0.fname ',1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm('defaults', 'EEG');
        spm_jobman('run', matlabbatch);
    %     spm_smooth(P0,fullfile(P,['sDelay_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s')  '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']),[3 3 3]);
        Q=fullfile(P,['sDelay_' NameEpileptogenicity '_' FileName spm_str_manip(Dinit.fname,'s')  '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']);
        M=spm_vol(Q);
        I=spm_read_vols(M);
        I(find(isnan(Delay)))=NaN;
        I=(max(Delay(:))./max(I(:)))*I;
        M = spm_write_vol(M,I);
    end
    if size(DD,1)>1
        Delay=NaN*zeros(53,63,46);
        Delay=NaN*zeros(V.dim(1),V.dim(2),V.dim(3));
        for i2=1:size(latency,2)
            Latency=mean(latency(:,i2));
            SPMFile=fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(Latency))], 'SPM.mat');
            load(SPMFile)
            % Load spmT map
            P1 = spm_vol(fullfile(P,['SPM_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(round(Latency))],'spmT_0001.nii'));
            Vb1 = spm_read_vols(P1);
            % Calcule threshold for statistical map
            df = [SPM.xCon(1).eidf SPM.xX.erdf];
            S    = SPM.xVol.S;                  %-search Volume {voxels}
            R    = SPM.xVol.R;                  %-search Volume {resels}
    %         u = spm_uc_FDR(0.001,df,'T',1,P1);
            u = spm_uc(ThDelay,df,'T',R,1,S);
            % Activated voxels
            Q1=find(Vb1>=u);
            Q2=find(isnan(Delay));
            Q3=intersect(Q2,Q1);
            Delay(Q3)=latency(i2);
        end
        P0=P1;
        P0.fname=fullfile(P,['Delay_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']);
        P0 = spm_write_vol(P0,Delay);
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {[P0.fname ',1']};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm('defaults', 'EEG');
        spm_jobman('run', matlabbatch);
    %     spm_smooth(P0,fullfile(P,['sDelay_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']),[3 3 3]);
        Q=fullfile(P,['sDelay_' NameEpileptogenicity '_Group_' FileName '_' num2str(min(FreqBand)) '_' num2str(max(FreqBand)) '_' num2str(round(mean(Horizon))) '_' num2str(1000*ThDelay) '.nii']);
        M=spm_vol(Q);
        I=spm_read_vols(M);
        I(find(isnan(Delay)))=NaN;
        I=(max(Delay(:))./max(I(:)))*I;
        M = spm_write_vol(M,I);
    end    
end