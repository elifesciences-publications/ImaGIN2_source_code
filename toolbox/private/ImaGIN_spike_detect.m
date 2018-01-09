function [SpikeRate,MinLoc,MaxLoc]=ImaGIN_spike_detect(Data,Time,FreqSpike,ThreshData,ThreshCC,FlagDetrend,Coarse,FlagContinuous)
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
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David

try
    FlagContinuous;
catch
    FlagContinuous=0;
end


%Data is a vector with artefact removed

%Downsample data
Time0=Time;
tmp=20*FreqSpike;
try
    Coarse;
catch
    Coarse=floor((1/(Time(2)-Time(1)))/tmp);
end
if isempty(Coarse)||Coarse==0
    Coarse=floor((1/(Time(2)-Time(1)))/tmp);
end    
Data=Data(1:Coarse:end);
Time=Time(1:Coarse:end);


TimeWindowWidth=1/FreqSpike;
TimeWindowWidth2=ceil(0.25*TimeWindowWidth/(Time(2)-Time(1)));
TimeWindowWidth=ceil(0.5*TimeWindowWidth/(Time(2)-Time(1)));
BigWin=20*TimeWindowWidth2;

%find local Maxima & Minima
MinLoc=zeros(size(Data));
MaxLoc=zeros(size(Data));
if ThreshData==0
    ThreshData=std(Data);
end
if ~FlagDetrend
    for i1=TimeWindowWidth2+1:length(Data)-TimeWindowWidth2-1
%         if Data(i1)==Data(i1-1)&Data(i1)==Data(i1+1)
%             if Data(i1)-Data(i1-2)<=0&Data(i1)-Data(i1+2)<=0&Data(i1)<-ThreshData&Data(i1)-Data(i1-TimeWindowWidth2)<=0&Data(i1)-Data(i1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if Data(i1)-Data(i1-2)>=0&Data(i1)-Data(i1+2)>=0&Data(i1)>ThreshData&Data(i1)-Data(i1-TimeWindowWidth2)>=0&Data(i1)-Data(i1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         else
%             if Data(i1)-Data(i1-1)<=0&Data(i1)-Data(i1+1)<=0&Data(i1)<-ThreshData&Data(i1)-Data(i1-TimeWindowWidth2)<=0&Data(i1)-Data(i1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if Data(i1)-Data(i1-1)>=0&Data(i1)-Data(i1+1)>=0&Data(i1)>ThreshData&Data(i1)-Data(i1-TimeWindowWidth2)>=0&Data(i1)-Data(i1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         end

%         if Data(i1)==Data(i1-1)&Data(i1)==Data(i1+1)
%             if Data(i1)-Data(i1-2)<=0&Data(i1)-Data(i1+2)<=0&Data(i1)-Data(i1-TimeWindowWidth2)<=0&Data(i1)-Data(i1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if Data(i1)-Data(i1-2)>=0&Data(i1)-Data(i1+2)>=0&Data(i1)-Data(i1-TimeWindowWidth2)>=0&Data(i1)-Data(i1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         else
%             if Data(i1)-Data(i1-1)<=0&Data(i1)-Data(i1+1)<=0&Data(i1)-Data(i1-TimeWindowWidth2)<=0&Data(i1)-Data(i1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if Data(i1)-Data(i1-1)>=0&Data(i1)-Data(i1+1)>=0&Data(i1)-Data(i1-TimeWindowWidth2)>=0&Data(i1)-Data(i1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         end

        if Data(i1)==Data(i1-1)&Data(i1)==Data(i1+1)
            if Data(i1)-Data(i1-2)<=0&Data(i1)-Data(i1+2)<=0
                MinLoc(i1)=1;
            end
            if Data(i1)-Data(i1-2)>=0&Data(i1)-Data(i1+2)>=0
                MaxLoc(i1)=1;
            end
        else
            if Data(i1)-Data(i1-1)<=0&Data(i1)-Data(i1+1)<=0
                MinLoc(i1)=1;
            end
            if Data(i1)-Data(i1-1)>=0&Data(i1)-Data(i1+1)>=0
                MaxLoc(i1)=1;
            end
        end


    end
else
%     for i1=BigWin+1:length(Data)-BigWin-1
%         win=i1-BigWin:i1+BigWin;
%         tmp=detrend(Data(win),2);
%         if tmp(BigWin+1)==tmp(BigWin+1-1)&tmp(BigWin+1)==tmp(BigWin+1+1)
%             if tmp(BigWin+1)-tmp(BigWin+1-2)<=0&tmp(BigWin+1)-tmp(BigWin+1+2)<=0&tmp(BigWin+1)<-ThreshData&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)<=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if tmp(BigWin+1)-tmp(BigWin+1-2)>=0&tmp(BigWin+1)-tmp(BigWin+1+2)>=0&tmp(BigWin+1)>ThreshData&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)>=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         else
%             if tmp(BigWin+1)-tmp(BigWin+1-1)<=0&tmp(BigWin+1)-tmp(BigWin+1+1)<=0&tmp(BigWin+1)<-ThreshData&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)<=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)<=0
%                 MinLoc(i1)=1;
%             end
%             if tmp(BigWin+1)-tmp(BigWin+1-1)>=0&tmp(BigWin+1)-tmp(BigWin+1+1)>=0&tmp(BigWin+1)>ThreshData&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)>=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)>=0
%                 MaxLoc(i1)=1;
%             end
%         end
%     end
    for i1=BigWin+1:length(Data)-BigWin-1
        win=i1-BigWin:i1+BigWin;
        tmp=detrend(Data(win),2);
        if tmp(BigWin+1)==tmp(BigWin+1-1)&tmp(BigWin+1)==tmp(BigWin+1+1)
            if tmp(BigWin+1)-tmp(BigWin+1-2)<=0&tmp(BigWin+1)-tmp(BigWin+1+2)<=0&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)<=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)<=0
                MinLoc(i1)=1;
            end
            if tmp(BigWin+1)-tmp(BigWin+1-2)>=0&tmp(BigWin+1)-tmp(BigWin+1+2)>=0&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)>=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)>=0
                MaxLoc(i1)=1;
            end
        else
            if tmp(BigWin+1)-tmp(BigWin+1-1)<=0&tmp(BigWin+1)-tmp(BigWin+1+1)<=0&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)<=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)<=0
                MinLoc(i1)=1;
            end
            if tmp(BigWin+1)-tmp(BigWin+1-1)>=0&tmp(BigWin+1)-tmp(BigWin+1+1)>=0&tmp(BigWin+1)-tmp(BigWin+1-TimeWindowWidth2)>=0&tmp(BigWin+1)-tmp(BigWin+1+TimeWindowWidth2)>=0
                MaxLoc(i1)=1;
            end
        end
    end
end
MinLoc=find(MinLoc);
MaxLoc=find(MaxLoc);

%Keep only large peaks
MinLoc2=zeros(size(MinLoc));
n=0;
for i1=1:length(MinLoc)
%     tmp=MaxLoc(find(max(find(MaxLoc<MinLoc(i1)))));
%     if Data(tmp)-Data(MinLoc(i1))>ThreshData
%         n=n+1;
%         MinLoc2(n)=MinLoc(i1);
%     else
%         tmp=MaxLoc(find(min(find(MaxLoc>MinLoc(i1)))));
%         if Data(tmp)-Data(MinLoc(i1))>ThreshData
%             n=n+1;
%             MinLoc2(n)=MinLoc(i1);
%         end
%     end
    tmp1=MaxLoc(max(find(MaxLoc<MinLoc(i1))));
    tmp2=MaxLoc(min(find(MaxLoc>MinLoc(i1))));
    if ~isempty(tmp1)&&~isempty(tmp2)
        if Data(tmp1)-Data(MinLoc(i1))>ThreshData && Data(tmp2)-Data(MinLoc(i1))>ThreshData
            n=n+1;
            MinLoc2(n)=MinLoc(i1);
        end
    end
end
MinLoc2=MinLoc2(1:n);
MaxLoc2=zeros(size(MaxLoc));
n=0;
for i1=1:length(MaxLoc)
%     tmp=MinLoc(find(max(find(MinLoc<MaxLoc(i1)))));
%     if -Data(tmp)+Data(MaxLoc(i1))>ThreshData
%         n=n+1;
%         MaxLoc2(n)=MaxLoc(i1);
%     else
%         tmp=MinLoc(find(min(find(MinLoc>MaxLoc(i1)))));
%         if -Data(tmp)+Data(MaxLoc(i1))>ThreshData
%             n=n+1;
%             MaxLoc2(n)=MaxLoc(i1);
%         end
%     end
    tmp1=MinLoc(max(find(MinLoc<MaxLoc(i1))));
    tmp2=MinLoc(min(find(MinLoc>MaxLoc(i1))));
    if ~isempty(tmp1)&&~isempty(tmp2)
        if -Data(tmp1)+Data(MaxLoc(i1))>ThreshData && -Data(tmp2)+Data(MaxLoc(i1))>ThreshData
            n=n+1;
            MaxLoc2(n)=MaxLoc(i1);
        end
    end
end
MaxLoc2=MaxLoc2(1:n);
MaxLoc=MaxLoc2;
MinLoc=MinLoc2;


%try to remove non spikes first by assuming a standard shape
MinLoc2=zeros(size(MinLoc));
n=0;
for i1=1:length(MinLoc)
    tmp=setdiff(find(abs(MinLoc-MinLoc(i1))<=TimeWindowWidth),i1);
    if isempty(tmp)|Data(MinLoc(i1))<=min(Data(MinLoc(tmp)))
        if i1>1
            if Data(MinLoc(i1))~=Data(MinLoc(i1-1))
                n=n+1;
                MinLoc2(n)=MinLoc(i1);
            end
        else
            n=n+1;
            MinLoc2(n)=MinLoc(i1);
        end
    end
end
MinLoc2=MinLoc2(1:n);
MaxLoc2=zeros(size(MaxLoc));
n=0;
for i1=1:length(MaxLoc)
    tmp=setdiff(find(abs(MaxLoc-MaxLoc(i1))<=TimeWindowWidth),i1);
    if isempty(tmp)|Data(MaxLoc(i1))>=max(Data(MaxLoc(tmp)))
        if i1>1
            if Data(MaxLoc(i1))~=Data(MaxLoc(i1-1))
                n=n+1;
                MaxLoc2(n)=MaxLoc(i1);
            end
        else
            n=n+1;
            MaxLoc2(n)=MaxLoc(i1);
        end
    end
end
MaxLoc2=MaxLoc2(1:n);

       figure
       plot(Time,Data)
       hold on
       line([Time(MinLoc);Time(MinLoc)],[min(Data(:));max(Data(:))]*ones(1,length(MinLoc)),'color','k')
       line([Time(MinLoc2);Time(MinLoc2)],[min(Data(:));max(Data(:))]*ones(1,length(MinLoc2)),'color','r')
       hold off

MaxLoc=MaxLoc2;
MinLoc=MinLoc2;
MaxLoc=MaxLoc(find(MaxLoc>TimeWindowWidth&MaxLoc<length(Data)-TimeWindowWidth));
MinLoc=MinLoc(find(MinLoc>TimeWindowWidth&MinLoc<length(Data)-TimeWindowWidth));

       
if FlagContinuous
    %try to fill gaps and make local corrections assuming continuous
    %stimulation
    ISI=mean(diff(MaxLoc));
    ISIerror=0.05*ISI;
    for i0=1:2  %Two passes
        MaxLocNew=[];
        for i1=2:length(MaxLoc)-1;
            if MaxLoc(i1+1)-MaxLoc(i1)>2*(ISI-ISIerror)
                MaxLocNew=[MaxLocNew MaxLoc(i1)+round(ISI)];
            end
            if MaxLoc(i1)-MaxLoc(i1-1)<ISI-ISIerror&MaxLoc(i1+1)-MaxLoc(i1)>ISI+ISIerror;
                MaxLoc(i1)=round((MaxLoc(i1-1)+MaxLoc(i1+1))/2);
            elseif MaxLoc(i1)-MaxLoc(i1-1)>ISI+ISIerror&MaxLoc(i1+1)-MaxLoc(i1)<ISI-ISIerror;
                MaxLoc(i1)=round((MaxLoc(i1-1)+MaxLoc(i1+1))/2);
            end
        end
        MaxLoc=sort([MaxLoc MaxLocNew]);
    end
    ISI=mean(diff(MinLoc));
    ISIerror=0.05*ISI;
    for i0=1:2  %Two passes
        MinLocNew=[];
        for i1=2:length(MinLoc)-1;
            if MinLoc(i1+1)-MinLoc(i1)>2*(ISI-ISIerror)
                MinLocNew=[MinLocNew MinLoc(i1)+round(ISI)];
            end
            if MinLoc(i1)-MinLoc(i1-1)<ISI-ISIerror&MinLoc(i1+1)-MinLoc(i1)>ISI+ISIerror;
                MinLoc(i1)=round((MinLoc(i1-1)+MinLoc(i1+1))/2);
            elseif MinLoc(i1)-MinLoc(i1-1)>ISI+ISIerror&MinLoc(i1+1)-MinLoc(i1)<ISI-ISIerror;
                MinLoc(i1)=round((MinLoc(i1-1)+MinLoc(i1+1))/2);
            end
        end
        MinLoc=sort([MinLoc MinLocNew]);
    end
end

%try to remove non spikes first by constructing a template
if ~isempty(ThreshCC)        
        TemplateMax=zeros(1,2*TimeWindowWidth+1);
        for i1=1:length(MaxLoc)
            TemplateMax=TemplateMax+Data(MaxLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]);
        end
        TemplateMax=TemplateMax/length(MaxLoc);
        TemplateMin=zeros(1,2*TimeWindowWidth+1);
        for i1=1:length(MinLoc)
            TemplateMin=TemplateMin+Data(MinLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]);
        end
        TemplateMin=TemplateMin/length(MinLoc);
        
        ok=1;
        n0=0;
        while ok
            n0=n0+1;
            C=zeros(1,length(MaxLoc));
            for i1=1:length(MaxLoc)
                cc=corrcoef(Data(MaxLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]),TemplateMax);
                C(i1)=cc(2);
            end
            ThreshCCInit=ThreshCC;
            if ThreshCCInit==0 && n0==1
                y=ImaGIN_kMeansCluster(C',2);
                c=zeros(1,2);
                c(1)=mean(y(find(y(:,end)==1)));
                c(2)=mean(y(find(y(:,end)==2)));
%                 ThreshCC=(c(1)+c(2))/2;
                ThreshCC=min(c);
            end
            N=find(C>=ThreshCC);
            if length(N)~=length(C)
                MaxLoc=MaxLoc(N);
                TemplateMax=zeros(1,2*TimeWindowWidth+1);
                for i1=1:length(MaxLoc)
                    TemplateMax=TemplateMax+Data(MaxLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]);
                end
                TemplateMax=TemplateMax/length(MaxLoc);
%                 plot(TemplateMax);drawnow;length(MaxLoc)
            else
                ok=0;
            end
        end
        ok=1;
        n0=0;
        while ok
            n0=n0+1;
            C=zeros(1,length(MinLoc));
            for i1=1:length(MinLoc)
                cc=corrcoef(Data(MinLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]),TemplateMin);
                C(i1)=cc(2);
            end
            ThreshCCInit=ThreshCC;
            if ThreshCCInit==0 && n0==1
                y=ImaGIN_kMeansCluster(C',2);
                c=zeros(1,2);
                c(1)=mean(y(find(y(:,end)==1)));
                c(2)=mean(y(find(y(:,end)==2)));
%                 ThreshCC=(c(1)+c(2))/2;
                ThreshCC=min(c);
            end
            N=find(C>=ThreshCC);
            if length(N)~=length(C)
                MinLoc=MinLoc(N);
                TemplateMin=zeros(1,2*TimeWindowWidth+1);
                for i1=1:length(MinLoc)
                    TemplateMin=TemplateMin+Data(MinLoc(i1)+[-TimeWindowWidth:TimeWindowWidth]);
                end
                TemplateMin=TemplateMin/length(MinLoc);
%                 plot(TemplateMin);drawnow;length(MinLoc)
            else
                ok=0;
            end
        end

%        figure
%        plot(Time,Data)
%        hold on
%        line([Time(MinLoc2);Time(MinLoc2)],[min(Data(:));max(Data(:))]*ones(1,length(MinLoc2)),'color','r')
%        line([Time(MinLoc);Time(MinLoc)],[min(Data(:));max(Data(:))]*ones(1,length(MinLoc)),'color','k')
%        hold off

end       

%Compute instantaneous spike rate
SpikeRateMax=zeros(size(Data));
for i1=1:length(MaxLoc)-1
    SpikeRateMax(MaxLoc(i1):MaxLoc(i1+1))=1/((Time(2)-Time(1))*(MaxLoc(i1+1)-MaxLoc(i1)));
end
SpikeRateMin=zeros(size(Data));
for i1=1:length(MinLoc)-1
    SpikeRateMin(MinLoc(i1):MinLoc(i1+1))=1/((Time(2)-Time(1))*(MinLoc(i1+1)-MinLoc(i1)));
end
%        SpikeRate=mean([SpikeRateMax;SpikeRateMin]);
SpikeRate=max([SpikeRateMax;SpikeRateMin],[],1);

if Coarse>1
    SpikeRate=interp1(Time,SpikeRate,Time0);
    MinLoc=MinLoc*Coarse;
    MaxLoc=MaxLoc*Coarse;
end
