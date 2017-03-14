function Dnew = ImaGIN_ReduceMontage(S)
% Perform a reduce montage from a set of defined electrodes


try
    Filename=S.Filename;
catch
    Filename = spm_select(Inf, '\.mat$', 'Select data file');
end

try
    FileOut=S.FileOut;
catch
    S.FileOut = Filename;
end

try
    Ref = S.Ref;
catch
    Ref = {};
end

if ~isempty(Ref)
    for i0=1:size(Filename,1)
        T=deblank(Filename(i0,:));
        D=spm_eeg_load(T);
        
        %     Cnames=chanlabels(D);
        Sensors=sensors(D,'EEG');
        Name=chanlabels(D);
        BadChannelsMono=badchannels(D);
        
        Index=zeros(1,length(Ref));
        for i2=1:length(Ref)
            for i1=1:length(Name)-1
                if strcmp(Ref{i2},Name{i1})
                    Index(i2)=i1;
                end
            end
        end
        Name=Name(Index);
        if ~isempty(Sensors)
            Sensors=Sensors(Index);
        end
                
        %New bad channels
        BadChannelsNew=[];
        for i1=1:length(BadChannelsMono)
            [a1,a2]=find(Index==BadChannelsMono(i1));
            BadChannelsNew=[BadChannelsBip a2'];
        end
        BadChannelsNew=unique(BadChannelsNew);
        
        
        
        Data=D(Index,:,:);
        
        %Save as a newfile
        try
            newname=FileOut;
            Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
        catch
            newname=['r' D.fname];
            Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
        end
        Dnew(:,:,:)=Data;
        Dnew=sensors(Dnew,'EEG',Sensors);
        Dnew=chanlabels(Dnew,1:size(Data,1),Name);
        Dnew=chantype(Dnew,1:size(Data,1),chantype(D,Index));
        if ~isempty(BadChannelsNew)
            Dnew = badchannels(Dnew,BadChannelsNew,1); % add badchannel index in meeg object
        end
        
        %     for i1=1:size(Data,1)
        %         Dnew=chantype(Dnew,i1,chantype(D,bipole(1,i1)));
        %     end
        save(Dnew);
        
        
    end
    
end
    
