function usSignal=ImaGIN_DownSample(Signal, fe, fs)

step1=round(fe/fs);
step2=ceil(fe/fs);

usSignal=zeros(step2,length(1:step1:length(Signal)));
for i=1:step2
%     s=Signal(1:step1:end);
%     i_sp1=floor(i/(fe/fs));
%     if length(s)>1
%         usSignal(i,:)=cat(1,s(i_sp1+1:end),s(1:i_sp1));
%     else
%         usSignal(i,:)=s;
%     end
%     Signal=cat(1,Signal(end), Signal(1:end-1));

    index=i+[0:step1:length(Signal)-i];
    usSignal(i,1:length(index))=Signal(index);
end

end