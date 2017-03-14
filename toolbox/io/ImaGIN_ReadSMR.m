function [data,Channel]=ReadSMR(DataFile,channelnumber,Coarse,ReadTime)

fid=fopen(DataFile);
data.data=[];
Channel=[];
data.name={};
n=0;
n2=0;

%In case signals are not sampled at the same frequency
f=0;
for i1=1:length(channelnumber)
    [d,header]=SONGetChannel(fid,channelnumber(i1));
    if header.sampleinterval>f
        f=header.sampleinterval;
    end    
end

for i1=1:length(channelnumber)
    [d,header]=ImaGIN_SONGetChannel(fid,channelnumber(i1),ReadTime);
    Coarse2=round(Coarse*f/header.sampleinterval);
    d=d';
    d=d(:);
    if ~isfield(d,'timings')
        if length(d)~=length(find(d==-1))&&isfield(header,'start')
            n=n+1;
            if n==1
                data.data=zeros(length(channelnumber),length(d(1:Coarse2:end)));
                data.time=header.start(1):header.sampleinterval*Coarse2:header.stop(end);
            end
            %         data.data(n,1:length(double(d(1:Coarse:end))))=double(d(1:Coarse:end));
            if Coarse2>1
try
    %                 d= ImaGIN_lowpassFilter(double(d),1/header.sampleinterval,1/(2*header.sampleinterval*Coarse2)-1);    %antialiasing,
                d= ImaGIN_bandpassFilter(double(d),1/header.sampleinterval,1,1/(2*header.sampleinterval*Coarse2)-1);    %antialiasing,
end
            end
            if i1==1
                data.data(n,:)=d(1:Coarse2:end);
            else
                if size(data.data,2)~=length(d(1:Coarse2:end))
                    data.data=data.data(:,1:length(d(1:Coarse2:end)));
                    data.time=data.time(1:length(d(1:Coarse2:end)));
                end
                data.data(n,:)=d(1:Coarse2:end);
            end
            Channel=[Channel i1];
            data.name{n}=header.title;
        end
    else    %spikes
        n2=n2+1;
        data.spike.timings{n2}=d.timings;
        data.spike.markers{n2}=d.markers;
        data.spike.adc{n2}=d.adc;
        data.spike.name{n2}=header.title;
        data.spike.sampleinterval=header.sampleinterval;
        data.spike.FileChannel(n2)=header.FileChannel;
    end
end
if n<size(data.data,1)
    tmp=data.data(1:n,:);
    data.data=tmp;
end
fclose(fid);