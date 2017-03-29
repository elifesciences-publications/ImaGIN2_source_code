function D = ImaGIN_causality_rescale(S)
% Time-rescale several causality files to normalise the length of events (e.g. seizures) before averaging

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
% Authors: Olivier David, 2008

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG causality rescale',0);

%files to rescale
try
    DD = S.D;
catch
    DD = spm_select(inf, '\.mat$', 'Select causality mat file(s) to rescale',[],pwd,'ca1_');
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
TimeTemplate2=TimeTemplate; %possible to resample EEG data differently 

for i2=1:length(DD)
    D=DD{i2};
    Time=D.ca.Time-S.Event(i2,1);
%     Time=Time/Time(find(abs(D.ca.Time-S.Event(i2,2))==min(abs(D.ca.Time-S.Event(i2,2)))));
    Time=Time/(S.Event(i2,2)-S.Event(i2,1));
    Time2=D.time-S.Event(i2,1);
    Time2=Time2/(S.Event(i2,2)-S.Event(i2,1));
    

    %Resample
    switch D.ca.Method
        case 'Autoregressive Models'

            %for surrogates
            Time_S=0:size(D.ca.S_S,3)-1;
            Time_S=min(TimeTemplate)+((max(TimeTemplate)-min(TimeTemplate))/max(Time_S))*Time_S;

            data=zeros(D.nchannels,length(TimeTemplate2));
            for i3=1:D.nchannels
                data(i3,:)=interp1(Time2,squeeze(D(i3,:)),TimeTemplate2);
            end
            s=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            s_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            dtf=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            dtf_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ddtf=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ddtf_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ffdtf=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ffdtf_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gddtf=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gddtf_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            pdc=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            pdc_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ffpdc=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            ffpdc_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gffpdc=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gffpdc_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gpdc=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            gpdc_S=zeros(D.nchannels,D.nchannels,length(TimeTemplate),length(D.ca.Freq));
            for i3=1:D.nchannels
                for i4=1:D.nchannels
                    for i5=1:length(D.ca.Freq)
                        s(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.S(i3,i4,:,i5)),TimeTemplate);
                        s_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.S_S(i3,i4,:,i5)),TimeTemplate);
                        dtf(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.DTF(i3,i4,:,i5)),TimeTemplate);
                        dtf_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.DTF_S(i3,i4,:,i5)),TimeTemplate);
                        ddtf(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.dDTF(i3,i4,:,i5)),TimeTemplate);
                        ddtf_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.dDTF_S(i3,i4,:,i5)),TimeTemplate);
                        ffdtf(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.ffDTF(i3,i4,:,i5)),TimeTemplate);
                        ffdtf_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.ffDTF_S(i3,i4,:,i5)),TimeTemplate);
                        gddtf(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.GdDTF(i3,i4,:,i5)),TimeTemplate);
                        gddtf_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.GdDTF_S(i3,i4,:,i5)),TimeTemplate);
                        pdc(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.PDC(i3,i4,:,i5)),TimeTemplate);
                        pdc_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.PDC_S(i3,i4,:,i5)),TimeTemplate);
                        ffpdc(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.ffPDC(i3,i4,:,i5)),TimeTemplate);
                        ffpdc_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.ffPDC_S(i3,i4,:,i5)),TimeTemplate);
                        gffpdc(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.GffPDC(i3,i4,:,i5)),TimeTemplate);
                        gffpdc_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.GffPDC_S(i3,i4,:,i5)),TimeTemplate);
                        gpdc(i3,i4,:,i5)=interp1(Time,squeeze(D.ca.GPDC(i3,i4,:,i5)),TimeTemplate);
                        gpdc_S(i3,i4,:,i5)=interp1(Time_S,squeeze(D.ca.GPDC_S(i3,i4,:,i5)),TimeTemplate);
                    end
                end
            end
            D.ca.DTF=dtf;
            D.ca.dDTF=ddtf;
            D.ca.PDC=pdc;
            D.ca.GPDC=gpdc;
            D.ca.DTF_S=dtf_S;
            D.ca.dDTF_S=ddtf_S;
            D.ca.PDC_S=pdc_S;
            D.ca.GPDC_S=gpdc_S;
            D1.ca.S=s;
            D1.ca.ffDTF=ffdtf;
            D1.ca.GdDTF=gddtf;
            D1.ca.ffPDC=ffpdc;
            D1.ca.GffPDC=gffpdc;
            D1.ca.S_S=s_S;
            D1.ca.ffDTF_S=ffdtf_S;
            D1.ca.GdDTF_S=gddtf_S;
            D1.ca.ffPDC=ffpdc;
            D1.ca.GffPDC_S=gffpdc_S;
    
            D.ca.p=interp1([1:length(D.ca.p)]/length(D.ca.p),D.ca.p,[1:length(TimeTemplate)]/length(TimeTemplate));
            D.ca.Order=interp1([1:length(D.ca.Order)]/length(D.ca.Order),D.ca.Order,[1:length(TimeTemplate)]/length(TimeTemplate));

        case 'Generalised Synchronisation'

            %for surrogates
            Time_S=0:size(D.ca.GSX_S,2)-1;
            Time_S=min(TimeTemplate)+((max(TimeTemplate)-min(TimeTemplate))/max(Time_S))*Time_S;

            data=zeros(D.nchannels,length(TimeTemplate2));
            for i3=1:D.nchannels
                data(i3,:)=interp1(Time2,squeeze(D(i3,:)),TimeTemplate2);
            end
            gsx=zeros(D.nchannels,length(TimeTemplate));
            gsx_S=zeros(D.nchannels,length(TimeTemplate));
            gsy=zeros(D.nchannels,length(TimeTemplate));
            gsy_S=zeros(D.nchannels,length(TimeTemplate));
            for i3=1:D.nchannels
                gsx(i3,:)=interp1(Time,D.ca.GSX(i3,:),TimeTemplate);
                gsx_S(i3,:)=interp1(Time_S,D.ca.GSX_S(i3,:),TimeTemplate);
                gsy(i3,:)=interp1(Time,D.ca.GSY(i3,:),TimeTemplate);
                gsy_S(i3,:)=interp1(Time_S,D.ca.GSY_S(i3,:),TimeTemplate);
            end
            D.ca.GSX=gsx;
            D.ca.GSX_S=gsx_S;
            D.ca.GSY=gsy;
            D.ca.GSY_S=gsy_S;
            D.ca.p=interp1([1:length(D.ca.p)]/length(D.ca.p),D.ca.p,[1:length(TimeTemplate)]/length(TimeTemplate));

    end

    %Write
    D=clone(D,['t' D.fnamedat],[
    D()=data;
    D.Nsamples=length(TimeTemplate2);
    D.ca.Time=TimeTemplate;
    D.time=TimeTemplate2;
    D=timeonset(D,find(abs(D.time)==min(abs(D.time))));
    D=fsample(D,1/(D.time(2)-D.time(1)));


    D.scale = spm_eeg_write(fpd, data, 2, D.datatype);

    save(D);
