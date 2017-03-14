function [data,Channel,evt]=ReadSMR(DataFile,channelnumber,Coarse)

%Works only if first electrode has one segment



fid=fopen(DataFile);
data.data=[];
Channel=[];
evt=[];
data.name={};
n=0;
n2=0;
n3=0;
n4=0;

%In case signals are not sampled at the same frequency
f=0;
Good=[];
Bad=[];
for i1=1:length(channelnumber)
    try
        [d,header]=SONGetChannel(fid,channelnumber(i1));
%         if ~isfield(header,'sampleinterval')
%             header.sampleinterval=1/50e3;
%         end
        if header.sampleinterval>1
            header.sampleinterval=1e-6*header.sampleinterval;
        end
        if header.sampleinterval>f
            f=header.sampleinterval;
        end
        if header.kind==3
            Bad=[Bad i1];
        else
            Good=[Good i1];
        end
    catch
        Bad=[Bad i1];
    end 
end

for i1=Good
    [d,header]=SONGetChannel(fid,channelnumber(i1));
%         if ~isfield(header,'sampleinterval')
%             header.sampleinterval=1/50e3;
%         end
    try
        if header.sampleinterval>1
            header.sampleinterval=1e-6*header.sampleinterval;
        end
    end
    try
        header.sampleinterval;
    catch
        header.sampleinterval=f;
    end

    
    %reconnect segments if interval is small
    try
        ok=1;
        while ok
            if length(header.start)==1
                ok=0;
            else
                interval=header.start(2:end)-header.stop(1:end-1);
                intervalsample=interval./header.sampleinterval;
                for i2=1:length(intervalsample)
                    if intervalsample(i2)<=10
                        ok=1;
                        Time1=[header.start(i2):header.sampleinterval:header.stop(i2)];
                        Time2=[header.start(i2+1):header.sampleinterval:header.stop(i2+1)];
                        TimeA=[header.start(i2):header.sampleinterval:header.stop(i2) header.start(i2+1):header.sampleinterval:header.stop(i2+1)];
                        TimeB=[header.start(i2):header.sampleinterval:header.stop(i2+1)];
                        Data1=double([d(i2,1:length(Time1)) d(i2+1,1:length(Time2))]);
                        Data2=interp1(TimeA,Data1,TimeB);
                        if length(header.start)<=i2+1
                            header.start=header.start(1:i2);
                        else
                            header.start=header.start([1:i2 i2+2:end]);
                        end
                        header.stop=header.stop([1:i2-1 i2+1:end ]);
                        d=d(setdiff(1:size(d,1),i2+1),:);
                        d(i2,1:length(Data2))=Data2;
                        if size(d,1)==1
                            d=d(:,1:length(TimeB));
                        end
                        break
                    else
                        ok=0;
                    end
                end
            end
        end
    end
                            
                
    
    
    %select longest segment
    try
        tmp=header.stop-header.start;
        if length(tmp)>1
            Segment=find(tmp==max(tmp));
            d=d(Segment,:);
        else
            Segment=1;
        end
    end



    Coarse2=round(Coarse*f/header.sampleinterval);
    d=d';
    d=d(:);
    if ~isfield(d,'timings')
        if length(d)~=length(find(d==-1))&&isfield(header,'start')
            n=n+1;
            if n==1
                data.data=zeros(length(channelnumber),length(d(1:Coarse2:end)));
%                 data.time=header.start(1):header.sampleinterval*Coarse2:header.stop(end);
                data.time=header.start(Segment):header.sampleinterval*Coarse2:header.stop(Segment);
                if size(data.data,2)>length(data.time)
                    data.time=0:length(d(1:Coarse2:end))-1;
                    data.time=header.sampleinterval*Coarse2*data.time+header.start(1);
                    data.time=header.start(Segment)+data.time;
                end
            end
            %         data.data(n,1:length(double(d(1:Coarse:end))))=double(d(1:Coarse:end));
            d=10*double(d)./(2^16-1);  %spike2 data are in int16
            if Coarse2>1
try
    %                 d= ImaGIN_lowpassFilter(double(d),1/header.sampleinterval,1/(2*header.sampleinterval*Coarse2)-1);    %antialiasing,
                d= ImaGIN_bandpassFilter(d,1/header.sampleinterval,1,1/(2*header.sampleinterval*Coarse2)-1);    %antialiasing,
end
            end
            if i1==1
                data.data(n,:)=d(1:Coarse2:end);
%                 data.time=header.start(Segment):header.sampleinterval*Coarse2:header.stop(Segment);
            else
                time=header.start(Segment):header.sampleinterval*Coarse2:header.stop(Segment);
                onset=min(find(abs(time-min(data.time))==min(abs(time-min(data.time)))));
                offset=min(find(abs(time-max(data.time))==min(abs(time-max(data.time)))));
                onset2=min(find(abs(data.time-min(time))==min(abs(data.time-min(time)))));
                offset2=min(find(abs(data.time-max(time))==min(abs(data.time-max(time)))));
                if onset==1 && offset==length(data.time)
                    tmp=d(onset:Coarse2:end);
                    if length(tmp)>length(data.time)
                        tmp=tmp(1:length(data.time));
                    end
                    data.data(n,onset:offset)=tmp;
                elseif onset>1
                    tmp=d(onset:Coarse2:end);
                    if length(tmp)>length(data.time)
                        tmp=tmp(1:length(data.time));
                    end
                    data.data(n,1:length(tmp))=tmp;
                else
                    tmp=d(1:Coarse2:end);
                    if length(tmp)>length(data.time)-onset2+1
                        tmp=tmp(1:length(data.time)-onset2+1);
                    end
                    data.data(n,onset2-1+[1:length(tmp)])=tmp;
                end
            end
            Channel=[Channel i1];
            data.name{n}=header.title;
        end
    else    %spikes
        if isfield(d,'adc')
            n2=n2+1;
            data.spike.timings{n2}=d.timings;
            data.spike.markers{n2}=d.markers;
            data.spike.adc{n2}=d.adc;
            data.spike.name{n2}=header.title;
            data.spike.sampleinterval=header.sampleinterval;
            data.spike.FileChannel(n2)=header.FileChannel;
        else
            n4=n4+1;
            for i2=1:length(d.timings)
                n3=n3+1;
                evt(n3).type  = header.title;
                evt(n3).time = d.timings(i2);
                evt(n3).value=n4;
            end
        end
    end
end

for i1=Bad
    [d,header]=SONGetChannel(fid,channelnumber(i1));
    if header.kind==4||header.kind==3||header.kind==2
        n4=n4+1;
        for i2=1:length(d)
            n3=n3+1;
            evt(n3).type  = header.title;
            evt(n3).time = d(i2);
            evt(n3).value=n4;
        end
    elseif header.kind==5
        n4=n4+1;
        for i2=1:size(d.timings,1)
            n3=n3+1;
            evt(n3).type  = header.title;
            evt(n3).time = d.timings(i2);
            evt(n3).value=n4;
        end
    end
end

    


if n<size(data.data,1)
    tmp=data.data(1:n,:);
    data.data=tmp;
end
fclose(fid);