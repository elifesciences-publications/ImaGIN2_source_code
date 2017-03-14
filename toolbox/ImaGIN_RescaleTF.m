function D = ImaGIN_RescaleTF(S)

%Time-rescale several TF files to normalise the length of events (e.g.
%seizures) before averaging
%Olivier David
%2008-03-14


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG time-frequency rescale',0);

%files to rescale
try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select TF EEG mat file(s) to rescale');
end

try
    clear tmp
    for i1=1:size(DD,1)
        tmp{i1} = spm_eeg_load(deblank(DD(i1,:)));
        Dname{i1}=DD(i1,:);
    end
    DD=tmp;
catch
    error(sprintf('Trouble reading file %s', DD));
end

%start and end of events
try
    S.Event;
catch
    for i2=1:length(DD)
        spm_input(sprintf('Onset/End for file %d',i2),1,'d');
        S.Event(i2,1) = spm_input('Onset [s]', '+1', 'r', '', 1);
        S.Event(i2,2) = spm_input('End [s]', '+1', 'r', '', 1);
    end
end
    

try
    S.Time
catch
    S.Time(1) = spm_input('Time onset (<0)', 1, 'r', '-0.2', 1);
    S.Time(2) = spm_input('Time end (>1)', '+1', 'r', '1.2', 1);
    S.Time(3) = spm_input('Time resolution (<1)', '+1', 'r', '0.001', 1);
end
TimeTemplate=[S.Time(1):S.Time(3):S.Time(2)];

for i2=1:length(DD)
    D=DD{i2};
    Time=D.tf.time-S.Event(i2,1);
%     Time=Time/Time(find(abs(D.tf.time-S.Event(i2,2))==min(abs(D.tf.time-S.Event(i2,2)))));
    Time=Time/(S.Event(i2,2)-S.Event(i2,1));

    %Resample
    data=zeros(D.nchannels,D.Nfrequencies,length(TimeTemplate)); 
    Dnew=clone(D, ['t' spm_str_manip(fnamedat(D),'t')], [D.nchannels D.Nfrequencies length(TimeTemplate) D.ntrials]);
    for i3=1:D.nchannels
        for i4=1:D.Nfrequencies
            data(i3,i4,:)=interp1(Time,squeeze(D(i3,i4,:)),TimeTemplate);
        end
    end
    Dnew(:,:,:)=data;
    Dnew.tf.Label=['Rescaled ' lower(D.tf.Label)];
    Dnew.tf.time=TimeTemplate;
    Dnew=fsample(Dnew,1/(Dnew.tf.time(2)-Dnew.tf.time(1)));
    save(Dnew);

end
