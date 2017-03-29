function D=ImaGIN_ArtefactCorrectionModel(S)
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
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
% 
% Author: Olivier David

try
    t=S.Filename;
catch
    t=spm_select(Inf, '\.mat$', 'Select data file');
end

D=spm_eeg_load(t);

try
    tms=S.tms;
catch
    tms=spm_input('Artefact total duration', '+1', 'i',1);
end

try
    rc1=S.rc1;
catch
    rc1=spm_input('Min RC value (sec)', 'r',10^(-6));
end

try
    rc2=S.rc2;
catch
    rc2=spm_input('Max RC value (sec)', 'r', 10^(-2));
end

try
    nb=S.nb;
catch
    nb=spm_input('Number of artefact shapes generated', 'r',300);
end

try
    mode=S.mode;
catch
    mode=spm_input('Artefact mode : "biphasic" or "monophasic"', 's','biphasic');
end

try
    DirScreenshots=S.DirScreenshots;
catch
    [dir,~,~]=fileparts(t);
end

try
    channels = S.channels;
    if isempty(channels)
        channels = 1:size(D,1);
    end
catch
    channels = spm_input('Select channel', 1, 'i');
end

if strcmp(mode,'biphasic')
    try
        break_duration=S.break;
    catch
        break_duration=spm_input('Break duration between the two pulses', 'r', 25*10^(-6));
    end
else break_duration=0;
end

try
    screenshot = S.screenshot;
catch
    spm_input('Save screenshots ?', '+1', 'Yes|No');
end

try
    estimateLength = S.artefactEstimation;
catch
    spm_input('Estimate the artifact length ?', '+1', 'Yes|No');
end

try
    adjustement = S.adjustement;
catch
    adjustement = 0;
end

try
    FileOut=S.FileOut;
catch
    FileOut=S.Filename;
end

saveb=0;     %save the rc value chosed for the correction

Dcorr = D(:,:);
Bad = badchannels(D);
SelectedChannels = setdiff(channels,Bad);
GoodChannels = setdiff(1:D.nchannels,Bad);

ev=events(D);
IndexEvent=[];
for i1=1:length(ev)
    if strcmp(ev(i1).type,'Stim')
        IndexEvent(end+1)=i1;
    end
end


if estimateLength
    Data=D(GoodChannels,:);
    d=[zeros(size(Data,1),1) abs(diff(Data,2,2)) zeros(size(Data,1),1)];
    d=ImaGIN_Normalisation(d,2);
    d=max(d);
    d=ImaGIN_Normalisation(d,2);
    
    template=0;
    for e=IndexEvent
        template=template+d(indsample(D,ev(e).time)+[-round(tms*2*D.fsample):round(tms*3*D.fsample)]);
    end
    
    template=template/length(IndexEvent);
    sup=find(template>max(template)/2);
    tms=((length(sup)-1)/D.fsample);
    tms_int=[0.8*tms tms 1.2*tms];  %estimated values tends to overestimate the time length
  
else tms_int=[tms 1.1*tms 1.25*tms];  %values specified on the recordings tends to underestimate the time length
end


windowCor=1.5*tms;   % time window of the artefact fitting : corresponds to at least 3 samples
if windowCor*D.fsample<3
    windowCor=3/D.fsample;
end

windowArt=5*windowCor;   % time window of the artefact (after each stimulation)
fe_ovsmpl=50000;    % initial sampling frequency
IntCor=round(D.fsample*windowCor);

if rc1<1/fe_ovsmpl
    rc1=1/fe_ovsmpl;
end


mArt=[];
for i0=1:length(tms_int)
    mArt=cat(2,mArt,ImaGIN_GenerateArtRange(rc1, rc2, nb, fe_ovsmpl, tms_int(i0), break_duration, windowArt, mode));  % generation of the different artifact shapes according to the minimal and maximal RC values
end

step1=round(fe_ovsmpl/fsample(D));
step2=ceil(fe_ovsmpl/fsample(D));
mmArt=zeros(step2,size(mArt,2),length(1:step1:size(mArt,1)));
for i=1:size(mArt,2)
    mmArt(:,i,:)=ImaGIN_DownSample(mArt(:,i),fe_ovsmpl,fsample(D));    % undersampling of every artifact shape
end

%correct by 3 samples (add zero) in case shift of artefact detection
mmArt2=cat(3,zeros(size(mmArt,1),size(mmArt,2)),mmArt);
mmArt2=mmArt2(:,:,1:size(mmArt,3));


if D.fsample>=2048
    mmArt3=cat(3,zeros(size(mmArt,1),size(mmArt,2)),mmArt2);
    mmArt3=mmArt3(:,:,1:size(mmArt,3));
    %     mmArt5=cat(3,zeros(size(mmArt,1),size(mmArt,2)),mmArt4);
    %     mmArt5=mmArt5(:,:,1:size(mmArt,3));
    mmArt=cat(1,mmArt,mmArt2,mmArt3);
else
    mmArt=cat(1,mmArt,mmArt2);
end



for c = SelectedChannels
    sig = [];
    ind1 = zeros(1,length(ev));
    ind2 = zeros(1,length(ev));
    countEv = 0;
    %selection of the part of the signal to be corrected
    for i = IndexEvent
        countEv = countEv+1;
        ind1(countEv) = indsample(D,ev(i).time)-adjustement;  %adjust the stimulation time start
        ind2(countEv) = indsample(D,ev(i).time+windowArt);
        intlength = ind2(1)-ind1(1);
        sig = cat(1,sig,D(c,ind1(countEv):ind1(countEv)+intlength-1));
    end
    m = min(size(sig,2),size(mmArt,3));
    [sig_corr,rc(c,:)] = ImaGIN_ArtefactFit(sig(:,1:m), mmArt(:,:,1:m), 'same', IntCor+1);
    
    startpoints=ind1(find(ind1));
    endpoints=ind1(find(ind1))+m-1;
    all=0;
    
    if screenshot
        [~, Name, ~] = fileparts(FileOut);
        if all
            hFig = figure;
            set(hFig, 'Position', [0 0 400 1100])
            plotNumb = 10;
            figNumb = 0;
            for i = 1:length(startpoints)
                if fix((i-1)/plotNumb) > fix((i-2)/plotNumb)
                    figNumb = figNumb +1;
                    print(hFig,fullfile(DirScreenshots,[Name '_c' int2str(c) '_' int2str(figNumb)]),'-dpng');
                    hFig = figure;
                    set(hFig, 'Position', [0 0 400 1100])
                end
                %generate screenshots
                subplot(plotNumb, 1, mod(i-1,plotNumb)+1)
                plot(sig_corr(i,:), 'g')
                hold on;plot(D(c,startpoints(i):endpoints(i)), 'k')
                %correct the signal
            end
            figNumb = figNumb +1;
            print(hFig,fullfile(DirScreenshots,[Name '_c' int2str(c) '_' int2str(figNumb)]),'-dpng');
        end
        %save screenshot of the mean
        hFig = figure;
        set(hFig, 'Position', [0 0 400 200])
        plot(mean(sig_corr,1),'g')
        hold on;plot(mean(sig(:,1:m)),'k')
        title([D.fname ', channel :' c ', sensor :' D.sensors('eeg').label(c)])
        print(hFig,fullfile(DirScreenshots,[Name '_mean_c' int2str(c)]),'-dpng');
        disp(c)
        close
        %correct the signal
        for i = 1:length(startpoints)
            Dcorr(c,startpoints(i):endpoints(i)) = sig_corr(i,:);
        end
    else
        %correct the signal
        for i = 1:length(startpoints)
            Dcorr(c,startpoints(i):endpoints(i)) = sig_corr(i,:);
        end
    end
        
end

if saveb == 1
    ind = (rc2-rc1)/nb;
    rc = rc*ind;
    
    Filename = ['rc_' D.fname '.txt'];
    fid = fopen(Filename,'w');
    fprintf(fid,'%f\n',rc);
    fclose(fid);
end

D = clone(D,FileOut, [D.nchannels D.nsamples D.ntrials]);
D(:,:,:) = Dcorr;
save(D);

S = [];
S.D = FileOut;
S.filter.band = 'low';
S.filter.type = 'butterworth';
S.filter.order = 5;
S.filter.dir = 'twopass';
S.filter.PHz = 90;
S.FileOut = FileOut;
D = ImaGIN_spm_eeg_filter(S);

end
