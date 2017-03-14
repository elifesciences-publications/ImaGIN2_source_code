function ImaGIN_BPFilter(P)

try
    t=P.Fname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end

F  = spm_figure('GetWin','Interactive');
figure(F);clf
try
    BP=P.BP;
catch
    BP = spm_input('Bandpass filter (e.g. 5 40)','1','r');
end

for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Data=D(:,:,:);
    for i1=1:D.nchannels
        for i2=1:D.ntrials 
            Data(i1,:,i2)=ImaGIN_bandpassFilter(squeeze(Data(i1,:,i2)),D.fsample,BP(1),BP(2));
        end
    end

    %Save as a newfile
    P = spm_str_manip(T, 'H');
    Prefix=[num2str(BP(1)) '-' num2str(BP(2))];
    Dnew = clone(D, [Prefix '_' fname(D)], [D.nchannels D.nsamples D.ntrials]);
    
    Dnew(:,:,:)=Data;

    save(Dnew);
end
