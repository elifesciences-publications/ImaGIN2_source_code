function ImaGIN_LowPassFilter(P)

try
    t=P.LFname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.05,'DesignMethod','butter');

for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Data=D(:,:,:);
    for i1=1:D.nchannels
        for i2=1:D.ntrials
            Data(i1,:,i2)=filtfilt(d1, squeeze(Data(i1,:,i2)));
        end
    end

    %Save as a newfile
    D=clone(D,['lpf_' D.fname], [D.nchannels D.nsamples D.ntrials]);
    D(:,:,:)=Data;
    save(D);
end