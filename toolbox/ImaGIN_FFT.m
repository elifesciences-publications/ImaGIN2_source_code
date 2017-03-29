function D = ImaGIN_FFT(S)
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
% Authors: Olivier David

try
    DD = S.D;
catch
    
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG FFT setup',0);
    
    DD = spm_select(inf, 'mat', 'Select EEG mat file');
end

for i1=1:size(DD,1)
    try
        D=ImaGIN_FFT_main(deblank(DD(i1,:)),S{i1});
    catch
        try
            D=ImaGIN_FFT_main(deblank(DD(i1,:)),S);
        catch
            D=ImaGIN_FFT_main(deblank(DD(i1,:)));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=ImaGIN_FFT_main(D,S)
    
P = spm_str_manip(D, 'H');

try
	D = spm_eeg_load(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

% Check for arguments in D
if ~isfield(D, 'fft')
	D.fft = [];
end

try
    ZeroCrossing=S.ZeroCrossing;
catch
    tmp=spm_input('Zero Crossing','+1','Yes|No');
    switch tmp
        case{'No'}
            ZeroCrossing=0;
        case{'Yes'}
            ZeroCrossing=1;
    end
end

% time=0:1/D.fsample:(D.nsamples-1)/D.fsample;
% %time=time-time(D.TimeZero);

Time1=time(D);

try
    Ratio = S.Ratio;
catch
    Ratio = spm_input('Cut-off ratio (median if empty)', '+1', 'r',[]);
end

try
    if isempty(S.TimeRange)
        TimeRange=[Time1(1) Time1(end)];
    else
        TimeRange=S.TimeRange;
    end
catch
    TimeRange = ...
        spm_input('Time interval [s] ([min max])', '+1', 'r',[Time1(1) Time1(end)]);
end

try
    FrequencyResolution=S.FrequencyResolution;
catch
    FrequencyResolution = ...
        spm_input('Frequency resolution (Hz)', '+1', 'r', '', [1, inf]);
end
Duration=1/FrequencyResolution;
Nsamples=ceil(D.fsample*Duration);
Time=[0:Nsamples-1]./D.fsample;
Freq = ImaGIN_Time2Freq(Time);
try
    if isempty(S.FrequencyRange)
        FrequencyRange=[Freq(2) max(Freq)];
    else
        FrequencyRange=S.FrequencyRange;
    end
catch
    FrequencyRange = ...
        spm_input('Frequency range [Hz] ([min max])', '+1', 'r',[Freq(2) max(Freq)]);
end
FreqIndex=find(Freq>=FrequencyRange(1)&Freq<=FrequencyRange(2));

F={};
for i0=1:size(TimeRange,1)
    TimeStart=min(find(Time1>=TimeRange(i0,1)));
    TimeEnd=max(find(Time1<=TimeRange(i0,2)));
    
    N=length(TimeStart:round(Nsamples/3):TimeEnd-Nsamples);
    if ZeroCrossing
        index=zeros(N,D.nchannels);
        n=0;
        for i1=TimeStart:round(Nsamples/3):TimeEnd-Nsamples
            n=n+1;
            index(n,:)=sum(abs(diff(sign(D(:,i1:i1+Nsamples-1))')));
        end
        index = min(index,[],2);
    else
        index=ones(N,1);
    end
    N=length(TimeStart:round(Nsamples/3):TimeEnd-Nsamples);
    f=zeros(N,D.nchannels,length(FreqIndex));
    n=0;
    n2=0;
    Win=hanning(Nsamples)';
    for i1=TimeStart:round(Nsamples/3):TimeEnd-Nsamples
        n=n+1;
        if index(n)>=max(index)/5
            n2=n2+1;
            tmp=abs(fft(D(:,i1:i1+Nsamples-1).*(ones(size(D,1),1)*Win),[],2));
            f(n2,:,:)=tmp(:,FreqIndex);
        end
    end
    f=f(1:n2,:,:);
    if isempty(Ratio)
        F{i0}=squeeze(median(f,1));
    else
%         ff=sort(f,1);
%         F{i0}=squeeze(ff(round(Ratio*size(ff,1)),:,:));
        [tmp1,tmp2]=sort(sum(sum(abs(f),2),3));
        F{i0}=squeeze(mean(f(tmp2(1:round(Ratio*length(tmp2))),:,:)));
    end
end

D.fft.Freq=Freq(FreqIndex);
D.fft.Data=F;
D.fft.Param.TimeRange=TimeRange;
D.fft.Param.FrequencyRange=FrequencyRange;
D.fft.Param.FrequencyResolution=FrequencyResolution;
D.fft.Param.n=n;

save(D);
