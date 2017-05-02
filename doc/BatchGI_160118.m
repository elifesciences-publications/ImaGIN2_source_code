clear all

Root='C:\Users\GIN11\Documents\Documents\PETHFO';
RootSEEG='C:\Users\GIN11\Documents\Documents\';
Root='/Users/odavid/Documents/Data/SEEG/GI';
RootSEEG='/Users/odavid/Documents/Data/SEEG';
RootAnat='/Users/odavid/Documents/Data/SEEG/Anat';
Root='/Volumes/Disk 2/Data/SEEG/GI';
RootAnat='/Users/odavid/Documents/Data/SEEG/Anat';
Root='/Users/admin/Documents/Data/SEEG/GI';
RootSEEG='/Users/admin/Documents/Data/SEEG';
RootAnat='/Users/admin/Documents/Data/SEEG/Anat';


%Nom du dossier du patient : Patient.Name
%Nom des fichiers de crise : Patient.File
%Nom des fichiers de baseline : Patient.Baseline
%Index des �lectrodes Bad : Patient.BadChannel
%Latences d'analyse : Patient.Latency
%Largeur de la fen�tre d'analyse : Patient.TimeConstant
%Lateralisation de l'implantation : Patient.Side ('L' ou 'R')
%Bande de frequence des HFOs : Patient.FreqBand

%Ici on rentre les infos des patients


Patient{1}.Name='Rot13Djo';
% Patient{1}.File{1}='Baseline 1';
% Patient{1}.File{2}='Baseline 2';
% Patient{1}.File{3}='G 1 1';
% Patient{1}.File{4}='G 1 2';
% Patient{1}.File{5}='G 1 3';
% Patient{1}.File{6}='G 1 4';
% Patient{1}.File{7}='G 1 5';
% Patient{1}.File{8}='G 1 6';
% Patient{1}.File{9}='G 1 7';
% Patient{1}.File{10}='G 1 8';
% Patient{1}.File{11}='G 2 1';
% Patient{1}.File{12}='G 2 2';
% Patient{1}.File{13}='G 2 3';
% Patient{1}.File{14}='NG 1 1';
% Patient{1}.File{15}='NG 1 2';
% Patient{1}.File{16}='NG 2 1';
% Patient{1}.File{17}='NG 2 2';
% Patient{1}.File{18}='NG 2 3';
Patient{1}.File{1}='G 1 1';
Patient{1}.Baseline{1}='Short_Baseline 1';
Patient{1}.File{2}='G 1 2';
Patient{1}.Baseline{2}='Short_Baseline 1';
Patient{1}.File{3}='G 1 3';
Patient{1}.Baseline{3}='Short_Baseline 1';
Patient{1}.File{4}='G 1 4';
Patient{1}.Baseline{4}='Short_Baseline 1';
Patient{1}.File{5}='G 1 5';
Patient{1}.Baseline{5}='Short_Baseline 1';
Patient{1}.File{6}='G 1 6';
Patient{1}.Baseline{6}='Short_Baseline 1';
Patient{1}.File{7}='G 1 7';
Patient{1}.Baseline{7}='Short_Baseline 1';
Patient{1}.File{8}='G 1 8';
Patient{1}.Baseline{8}='Short_Baseline 1';
Patient{1}.File{9}='G 2 1';
Patient{1}.Baseline{9}='Short_Baseline 2';
Patient{1}.File{10}='G 2 2';
Patient{1}.Baseline{10}='Short_Baseline 2';
Patient{1}.File{11}='G 2 3';
Patient{1}.Baseline{11}='Short_Baseline 2';
Patient{1}.File{12}='NG 1 1';
Patient{1}.Baseline{12}='Short_Baseline 2';
Patient{1}.File{13}='NG 1 2';
Patient{1}.Baseline{13}='Short_Baseline 2';
Patient{1}.File{14}='NG 2 1';
Patient{1}.Baseline{14}='Short_Baseline 2';
Patient{1}.File{15}='NG 2 2';
Patient{1}.Baseline{15}='Short_Baseline 2';
Patient{1}.File{16}='NG 2 3';
Patient{1}.Baseline{16}='Short_Baseline 2';
Patient{1}.BadChannel=[46];
Patient{1}.FreqBand=[80 120];
Patient{1}.Latency=[0];
Patient{1}.TimeConstant=1;
Patient{1}.Side='L';
Patient{1}.CreateTemplate='No';
Patient{1}.Template='Rot13DjoMNI';
Patient{1}.sMRI=fullfile(RootAnat,Patient{1}.Name,'MRI','wT1.img');

% Patient{2}.Name='Rot13Poil';
% % Patient{2}.File{1}='Seizure1';
% % Patient{2}.Baseline{1}='Baseline_Seizure1';
% Patient{2}.File{1}='Seizure2';
% Patient{2}.Baseline{1}='Baseline_Seizure2';
% Patient{2}.File{2}='Seizure3';
% Patient{2}.Baseline{2}='Baseline_Seizure3';
% Patient{2}.BadChannel=[];
% Patient{2}.FreqBand=[60 100];
% Patient{2}.Latency=[0];
% Patient{2}.TimeConstant=1;
% Patient{2}.Template='Rot13PoilMNI';
% Patient{2}.sMRI=fullfile(RootAnat,Patient{2}.Name,'MRI','wBrainPre.nii');

Patient{2}.Name='Rot13Poil';
% Patient{2}.File{1}='HFO';
Patient{2}.File{1}='Onset_4-4_016HFO';
Patient{2}.Baseline{1}='Baseline';
for i1=1:16
    if i1<10
        Patient{2}.File{i1}=['Onset_4-4_00' num2str(i1) 'HFO'];
    else
        Patient{2}.File{i1}=['Onset_4-4_0' num2str(i1) 'HFO'];
    end
    Patient{2}.Baseline{i1}='Baseline';
end
Patient{2}.BadChannel=[65 66];
Patient{2}.FreqBand=[80 110];
Patient{2}.Latency=[0];
Patient{2}.TimeConstant=0.6;
Patient{2}.Template='Rot13PoilMNI';
Patient{2}.sMRI=fullfile(RootAnat,Patient{2}.Name,'MRI','wBrainPre.nii');

Patient{3}.Name='Rot13Car';
% Patient{3}.File{1}='Spasm';
% Patient{3}.Baseline{1}='Baseline';
for i1=1:16
    if i1<10
        Patient{3}.File{i1}=['Onset_4-4_00' num2str(i1) 'Spasm'];
    else
        Patient{3}.File{i1}=['Onset_4-4_0' num2str(i1) 'Spasm'];
    end
    Patient{3}.Baseline{i1}='Baseline';
end
Patient{3}.BadChannel=[54];
Patient{3}.FreqBand=[60 120];
Patient{3}.Latency=[0];
Patient{3}.TimeConstant=0.8;
Patient{3}.Template='Rot13CarMNI';
Patient{3}.sMRI=fullfile(RootAnat,Patient{3}.Name,'MRI','wBrainPre.nii');

Patient{4}.Name='Rot15DeW';
Patient{4}.File{1}='Seizure1';
Patient{4}.Baseline{1}='Baseline_Seizure1';
Patient{4}.BadChannel=[34 35];
Patient{4}.FreqBand=[60 120];
Patient{4}.Latency=[0];
Patient{4}.TimeConstant=4;
Patient{4}.Template='Rot15DeWMNI';
Patient{4}.sMRI=fullfile(RootAnat,Patient{4}.Name,'MRI','wBrainPre.nii');

Patient{5}.Name='Rot14Gud';
% Patient{5}.File{1}='Spasm';
% Patient{5}.Baseline{1}='Baseline';
for i1=1:5
    if i1<10
        Patient{5}.File{i1}=['Onset_4-4_00' num2str(i1) 'Spasm'];
    else
        Patient{5}.File{i1}=['Onset_4-4_0' num2str(i1) 'Spasm'];
    end
    Patient{5}.Baseline{i1}='Baseline';
end
Patient{5}.BadChannel=61;
Patient{5}.FreqBand=[80 120];
Patient{5}.Latency=[0];
Patient{5}.TimeConstant=0.8;
Patient{5}.Template='Rot14GudMNI';
Patient{5}.sMRI=fullfile(RootAnat,Patient{5}.Name,'MRI','wBrainPre.nii');

Patient{6}.Name='Rot13Hai';
% Patient{6}.File{1}='Spasm';
% Patient{6}.Baseline{1}='Baseline';
for i1=1:15
    if i1<10
        Patient{6}.File{i1}=['Onset_4-4_00' num2str(i1) 'Spasm'];
    else
        Patient{6}.File{i1}=['Onset_4-4_0' num2str(i1) 'Spasm'];
    end
    Patient{6}.Baseline{i1}='Baseline';
end
Patient{6}.BadChannel=[21 32 65];
Patient{6}.FreqBand=[80 120];
Patient{6}.Latency=[0];
Patient{6}.TimeConstant=0.8;
Patient{6}.Template='Rot13HaiMNI';
Patient{6}.sMRI=fullfile(RootAnat,Patient{6}.Name,'MRI','wBrainPre.nii');


I=7;
Patient{I}.Name='0007ROT';
% Patient{6}.File{1}='Spasm';
% Patient{6}.Baseline{1}='Baseline';
%Exclude File 1 because not the same montage and too little overlap
for i1=2:4
    if i1<10
        Patient{I}.File{i1-1}=[Patient{I}.Name '_' num2str(i1)];
    else
        Patient{I}.File{i1-1}=['0007ROT_0' num2str(i1)];
    end
    Patient{7}.Baseline{i1-1}=['Baseline_' Patient{I}.File{i1-1}];
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[70 110];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=0.8;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=8;
Patient{I}.Name='0018ROT';
for i1=1:2
    if i1<10
        Patient{I}.File{i1}=[Patient{I}.Name '_' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_' num2str(i1)];
    else
        Patient{I}.File{i1}=[Patient{I}.Name '_0' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_0' num2str(i1)];
    end
end
Patient{I}.BadChannel=[16];
Patient{I}.FreqBand=[80 120];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=0.8;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=9;
N=[4 2 6];
Patient{I}.Name='RotErves';
for i1=1:3
    if i1<10
        Patient{I}.File{i1}=[Patient{I}.Name '_' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_' num2str(i1)];
    else
        Patient{I}.File{i1}=[Patient{I}.Name '_0' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_0' num2str(i1)];
    end
end
n=0;
for i1=1:3
    for i2=1:N(i1)
        n=n+1;
        Patient{I}.File{n}=['Spasm_5-7_00' num2str(i2) Patient{I}.Name '_' num2str(i1)];
        Patient{I}.Baseline{n}=['Spasm_6--1_00' num2str(i2) Patient{I}.Name '_' num2str(i1)];
    end
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[90 120];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=1;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=10;
Patient{I}.Name='0060ROT';
for i1=1:3
    if i1<10
        Patient{I}.File{i1}=[Patient{I}.Name '_' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_' num2str(i1)];
    else
        Patient{I}.File{i1}=[Patient{I}.Name '_0' num2str(i1)];
        Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_0' num2str(i1)];
    end
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[100 140];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=1;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=11;
N=[2 1 6 2];
Patient{I}.Name='0061ROT';
% for i1=1:4
%     if i1<10
%         Patient{I}.File{i1}=[Patient{I}.Name '_' num2str(i1)];
%         Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_' num2str(i1)];
%     else
%         Patient{I}.File{i1}=[Patient{I}.Name '_0' num2str(i1)];
%         Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.Name '_0' num2str(i1)];
%     end
% end
n=0;
for i1=1:4
    for i2=1:N(i1)
        n=n+1;
        Patient{I}.File{n}=['Spasm_5-7_00' num2str(i2) Patient{I}.Name '_' num2str(i1)];
        Patient{I}.Baseline{n}=['Spasm_6--1_00' num2str(i2) Patient{I}.Name '_' num2str(i1)];
        if i1==4
            Patient{I}.Baseline{n}=['Spasm_10--5_00' num2str(i2) Patient{I}.Name '_' num2str(i1)];
        end            
    end
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[100 120];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=1;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=12;
Patient{I}.Name='0001NAN';
Patient{I}.File{1}='0001NAN_Seizure_1';
Patient{I}.File{2}='Onset_140-500_0010001NAN_Seizure_2';
Patient{I}.File{3}='Onset_140-500_0020001NAN_Seizure_2';
Patient{I}.File{4}='Onset_140-500_0030001NAN_Seizure_2';
for i1=1:length(Patient{I}.File)
	Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.File{i1}];
end
Patient{I}.BadChannel=[164:176];
Patient{I}.FreqBand=[100 200];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=8;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');

I=13;
Patient{I}.Name='0100ROT';
Patient{I}.File{1}='0100ROT_1';
Patient{I}.File{2}='0100ROT_2';
for i1=1:length(Patient{I}.File)
	Patient{I}.Baseline{i1}=['Baseline_' Patient{I}.File{i1}];
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[120 180];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=2;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');
Patient{I}.Pre='Onset';

I=14;
Patient{I}.Name='0100ROT';
Patient{I}.File{1}='SP_20-30_001_0100ROT_1';
Patient{I}.File{2}='SP_20-30_002_0100ROT_1';
Patient{I}.File{3}='SP_20-30_003_0100ROT_1';
Patient{I}.File{4}='SP_20-30_004_0100ROT_1';
Patient{I}.File{5}='SP_20-30_005_0100ROT_1';
Patient{I}.File{6}='SP_20-30_006_0100ROT_1';
Patient{I}.File{7}='SP_20-30_001_0100ROT_2';
Patient{I}.File{8}='SP_20-30_002_0100ROT_2';
Patient{I}.File{9}='SP_20-30_003_0100ROT_2';
for i1=1:6%length(Patient{I}.File)
	Patient{I}.Baseline{i1}=['Baseline_0100ROT_1'];
end
for i1=7:9%length(Patient{I}.File)
	Patient{I}.Baseline{i1}=['Baseline_0100ROT_2'];
end
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[120 180];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=2;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');
Patient{I}.Pre='SP';

I=15;
Patient{I}.Name='0117ROT_2016-10-2516';
Patient{I}.File{1}='Seizure1';
Patient{I}.File{2}='Seizure2';
Patient{I}.File{1}='mod OM 1-4_60-60_001_Seizure1';
Patient{I}.Baseline{1}='mod OM 1-4_40--20_001_Seizure1';
Patient{I}.File{2}='crise_60-60_001_Seizure2';
Patient{I}.Baseline{2}='crise_50--30_001_Seizure2';
Patient{I}.File{3}='crise_60-60_002_Seizure2';
Patient{I}.Baseline{3}='crise_50--30_002_Seizure2';
Patient{I}.File{4}='crise_60-60_003_Seizure2';
Patient{I}.Baseline{4}='crise_50--30_003_Seizure2';
Patient{I}.File{5}='crise_60-60_004_Seizure2';
Patient{I}.Baseline{5}='crise_50--30_004_Seizure2';
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[140 220];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=4;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','Pre','w0117ROT_25102016.nii');
Patient{I}.Pre='';

I=16;
Patient{I}.Name='0011HUH';
Patient{I}.File{1}='SZ3';
Patient{I}.File{2}='SZ4';
Patient{I}.File{3}='SZ5';
Patient{I}.File{4}='SZ12';
Patient{I}.Baseline{1}='Baseline_SZ3';
Patient{I}.Baseline{2}='Baseline_SZ4';
Patient{I}.Baseline{3}='Baseline_SZ5';
Patient{I}.Baseline{4}='Baseline_SZ12';
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[60 120];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=4;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','wBrainPre.nii');
Patient{I}.Pre='';

I=17;
Patient{I}.Name='0011HUH_25112016';
Patient{I}.File{1}='SZ01';
Patient{I}.File{2}='SZ02';
Patient{I}.File{3}='SZ0304';
Patient{I}.File{4}='SZ05';
Patient{I}.Baseline{1}='Baseline_SZ01';
Patient{I}.Baseline{2}='Baseline_SZ01';
Patient{I}.Baseline{3}='Baseline_SZ0304';
Patient{I}.Baseline{4}='Baseline_SZ05';
Patient{I}.BadChannel=[];
Patient{I}.FreqBand=[100 180];
Patient{I}.FreqBand=[210 230];
Patient{I}.Latency=[0];
Patient{I}.TimeConstant=2;
Patient{I}.sMRI=fullfile(RootAnat,Patient{I}.Name,'MRI','Pre','w0011HUH_25112016.nii');
Patient{I}.Pre='';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%En-dessous, on ne touche plus

ThDelay=0.05;

for i0=17%[3 5 6]%2:length(Patient)
    
    %Convert TRC to SPM
%     [files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name,'Crises'),'.TRC');
    for i1=1:length(Patient{i0}.File)
        
        %Convert TRC to SPM
        try
            clear S
            S.dataset=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.TRC']);
            S.Atlas='Human';
            S.channel=[];
            S.coarse=1;
            S.FileOut=[Patient{i0}.File{i1}];
            S.loadevents='yes';
    %         S.NeventType=0;
            D = ImaGIN_spm_eeg_converteeg2mat(S);
        catch
            clear S
            S.dataset=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.edf']);
            S.Atlas='Human';
            S.channel=[];
            S.coarse=1;
            S.FileOut=[Patient{i0}.File{i1}];
            S.SEEG='Yes';
            S.SizeMax=1e10;
            D = ImaGIN_spm_eeg_converteeg2mat(S);
        end
                    
        %Add electrode position
        clear S
        S.Fname=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']);
        S.filenamePos=fullfile(RootAnat,Patient{i0}.Name,'Implantation/Electrodes_Pos_MNI.txt');
        S.filenameName=fullfile(RootAnat,Patient{i0}.Name,'Implantation/Electrodes_Name.txt');
        S.FileOut=S.Fname;
        D = ImaGIN_Electrode(S);
        
        %Longitudinal bipolar montage
        clear S
        S.Fname=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']);
        S.SaveFile=Patient{i0}.File{i1};
        S.FileOut=S.Fname;
        D = ImaGIN_BipolarMontage(S);

    end

    for i1=1:length(Patient{i0}.Baseline)
        
        %Convert TRC to SPM
        clear S
        S.dataset=fullfile(Root,Patient{i0}.Name,[Patient{i0}.Baseline{i1} '.TRC']);
        S.Atlas='Human';
        S.channel=[];
        S.coarse=1;
        S.SaveFile=Patient{i0}.Baseline{i1};
        S.loadevents='yes';
%         S.NeventType=0;
        D = ImaGIN_spm_eeg_converteeg2mat(S);

        %Add electrode position
        clear S
        S.Fname=fullfile(Root,Patient{i0}.Name,[Patient{i0}.Baseline{i1} '.mat']);
        S.filenamePos=fullfile(RootAnat,Patient{i0}.Name,'Implantation/Electrodes_Pos_MNI.txt');
        S.filenameName=fullfile(RootAnat,Patient{i0}.Name,'Implantation/Electrodes_Name.txt');
        D = ImaGIN_Electrode(S);
        
        %Longitudinal bipolar montage
        clear S
        S.Fname=fullfile(Root,Patient{i0}.Name,[Patient{i0}.Baseline{i1} '.mat']);
        S.SaveFile=Patient{i0}.Baseline{i1};
        S.FileOut=S.Fname;
        D = ImaGIN_BipolarMontage(S);
        
    end
    
    
end

%add here marker of seizure onset


for i0=12%[3 5 6]%2:length(Patient)
    
    
    %Set bad channels
    for i1=1:length(Patient{i0}.File)
        D=spm_eeg_load(fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']));
        D=badchannels(D,Patient{i0}.BadChannel,1);
        save(D);
    end
    for i1=1:length(Patient{i0}.Baseline)
        D=spm_eeg_load(fullfile(Root,Patient{i0}.Name,[Patient{i0}.Baseline{i1} '.mat']));
        D=badchannels(D,Patient{i0}.BadChannel,1);
        save(D);
    end
    
end

% %Compute wavelet
% for i0=2:length(Patient)
%     cd(fullfile(Root,Patient{i0}.Name))
%     TimeResolutionTF=0.05;
%     TimeResolution=0.1;
%     for i1=1:length(Patient{i0}.File)
%         clear SS
%         SS.D=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']);
%         D=spm_eeg_load(SS.D);
%         SS.Synchro='No';
%         SS.Pre='';
%         SS.Method='Morlet wavelet';
%         SS.frequencies=10:3:230;
%         SS.FactMod=0;
%         SS.Mfactor=20;
%         SS.Width=0;
%         SS.TimeWindow=-10:.2:10;
%         SS.TimeWindowWidth=SS.TimeWindow(2)-SS.TimeWindow(1);
%         SS.Coarse=0;
%         SS.channels=1:D.nchannels;
%         SS.TimeResolution=TimeResolutionTF;
%         ImaGIN_spm_eeg_tf(SS);
%         clear SS2
%         SS2.D=fullfile(D.path,['w1_' SS.Pre '_' D.fname]);
%         SS2.B=[-10 -1];
%         ImaGIN_NormaliseTF(SS2);
%     end
%     
%     
%     [files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name),'^nw1.*\.mat$');
%     if size(files,1)>1
%             cd(fullfile(Root,Patient{i0}.Name,'Crises'))
%         clear S
%     %     for i1=1:size(files,1)
%     %         S.D{i1}=deblank(files(i1,:));
%     %     end
%         S.D=files;
%         S.Method='Mean';
%         S.NewName='Mean';
%         D = ImaGIN_AverageTF(S);
%     end
% end
    
    
    
% Compute multi-taper
for i0=17%5:length(Patient)
    cd(fullfile(Root,Patient{i0}.Name))
    maxTime=10;
%     maxTime=2;
    minTime=-10;
    for i1=1:length(Patient{i0}.File)
        clear SS
        SS.D=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']);
        D=spm_eeg_load(SS.D);
        minTime=max([minTime min(time(D))+0.5]);
        maxTime=min([maxTime max(time(D))-0.5]);
    end
    for i1=1:length(Patient{i0}.File)
        clear SS
        SS.D=fullfile(Root,Patient{i0}.Name,[Patient{i0}.File{i1} '.mat']);
        D=spm_eeg_load(SS.D);
        SS.Pre='';
        SS.Method='Multitaper';
        SS.Taper='Hanning';
        SS.TimeResolution=0.1;
        SS.frequencies=10:3:230;
        SS.FactMod=10;
        SS.TimeWindow=[minTime maxTime];
        SS.TimeWindowWidth=1;
        SS.channels=1:D.nchannels;
        SS.NSegments=1;
        ImaGIN_spm_eeg_tf(SS);
        clear SS2
        SS2.D=fullfile(D.path,['m1_' SS.Pre '_' D.fname]);
        SS2.B=[-10 -1];
        ImaGIN_NormaliseTF(SS2);
    end
    
    
    [files,dirs] = spm_select('List',fullfile(Root,Patient{i0}.Name),'^nm1.*\.mat$');
    if size(files,1)>1
            cd(fullfile(Root,Patient{i0}.Name))
        clear S
    %     for i1=1:size(files,1)
    %         S.D{i1}=deblank(files(i1,:));
    %     end
        S.D=files;
        S.Method='Mean';
        S.NewName='Mean';
        D = ImaGIN_AverageTF(S);
    end
end


for i0=17:length(Patient)
    clear S
    for i1=1:length(Patient{i0}.File)
        if i1==1
            S.D=fullfile(Root,[Patient{i0}.Name],[Patient{i0}.File{i1} '.mat']);
            S.B=fullfile(Root,[Patient{i0}.Name],[Patient{i0}.Baseline{i1} '.mat']);
        else
            S.D=char(S.D,fullfile(Root,[Patient{i0}.Name],[Patient{i0}.File{i1} '.mat']));
            S.B=char(S.B,fullfile(Root,[Patient{i0}.Name],[Patient{i0}.Baseline{i1} '.mat']));
        end
    end
    S.TimeWindow=[0:0.01:Patient{i0}.TimeConstant+1+max(Patient{i0}.Latency)];
    S.FreqBand=Patient{i0}.FreqBand;
    S.HorizonT=Patient{i0}.TimeConstant;
    S.BadChannel=Patient{i0}.BadChannel;
    try
        S.FileName=Patient{i0}.Pre;
    catch
        S.FileName='';
    end
    S.Latency=Patient{i0}.Latency;
    S.TimeResolution=0.1;
    S.ThDelay=ThDelay;
    S.Atlas='Human';
    S.AR=0;
    S.Latency=0;
    S.sMRI=Patient{i0}.sMRI;
    S.CorticalMesh=1;
    ImaGIN_Epileptogenicity(S);
end