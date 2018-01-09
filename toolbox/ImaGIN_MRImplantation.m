function ImaGIN_MRImplantation(S)
% Define SEEG implantation
%
% DOCUMENTATION:
%    Creates 2 pts files (in original coordinates and in MNI space) using the .txt file 
%    created by looking at original structural patient's MRI. 
%    Needs the sn.mat created when normalising patient's MRI.

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
% Authors: Olivier David, 2007


try
    FileName=S.FileName;
catch
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','SEEG/MRI implantation setup',0);
    FileName = spm_select(1, '\.txt$', 'Select .txt file of electrodes');
end


[pathFile,nameFile,~]=fileparts(FileName);

try
    DirOut=S.DirOut;
catch
    DirOut=pathFile;
end

try
    Deformation=S.Deformation;
catch
    Deformation = spm_select(1, 'y_', 'Select forward deformation field y_');
end

try
    FileSource=S.FileSource;
catch
    tmp1=spm_str_manip(Deformation,'t');
    tmp2=spm_str_manip(Deformation,'h');
    FileSource=fullfile(tmp2,tmp1(3:end));
    if ~exist(FileSource)
        if ~exist([FileSource(1:end-3) 'img'])
            FileSource = spm_select(1, 'image', 'Select source image');
        else
            FileSource=[FileSource(1:end-3) 'img'];
        end
    end
end

NameElec={};
PosElec=[];
NName=0;
NPos=0;
d=3.5;
fid=fopen(FileName);
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    else
        if length(tline)>6
            NPos=NPos+1;
            PosElec(NPos,:)=str2num(tline);
        else
            NName=NName+1;
            NameElec{NName}=tline;
        end
    end
end
fclose(fid);

NElectrode=0;
NameElectrode={};
PosElectrode=[];
for i1=1:NName
    Pos=PosElec((i1-1)*2+[1:2],:);
    X = Pos(1,1):(Pos(2,1)-Pos(1,1))/100:Pos(2,1);
    Y = Pos(1,2):(Pos(2,2)-Pos(1,2))/100:Pos(2,2);
    Z = Pos(1,3):(Pos(2,3)-Pos(1,3))/100:Pos(2,3);
    if isempty(X)
        X=ones(1,101)*Pos(1,1);
    end
    if isempty(Y)
        Y=ones(1,101)*Pos(1,2);
    end
    if isempty(Z)
        Z=ones(1,101)*Pos(1,3);
    end
    D=sqrt((X-Pos(1,1)).^2+(Y-Pos(1,2)).^2+(Z-Pos(1,3)).^2);
    NElec=ceil(max(D)/d);
    for i2=1:NElec
        NElectrode=NElectrode+1;
        NameElectrode{NElectrode}=[NameElec{i1} num2str(i2)];
        tmp=max(find(D<=(i2-1)*d));
        PosElectrode(NElectrode,1:3)=[X(tmp) Y(tmp) Z(tmp)];
    end
end

%longitudinal bipolar montage
NElectrodeBip=0;
PosElectrodeBip=[];
NameElectrodeBip={};
bipole=[];
for i1=1:length(NameElectrode)-1
    tmp1=NameElectrode{i1}(1:min([findstr(NameElectrode{i1},'0') findstr(NameElectrode{i1},'1') findstr(NameElectrode{i1},'2') findstr(NameElectrode{i1},'3') findstr(NameElectrode{i1},'4') findstr(NameElectrode{i1},'5') findstr(NameElectrode{i1},'6') findstr(NameElectrode{i1},'7') findstr(NameElectrode{i1},'8') findstr(NameElectrode{i1},'9')])-1);
    tmp2=NameElectrode{i1+1}(1:min([findstr(NameElectrode{i1+1},'0') findstr(NameElectrode{i1+1},'1') findstr(NameElectrode{i1+1},'2') findstr(NameElectrode{i1+1},'3') findstr(NameElectrode{i1+1},'4') findstr(NameElectrode{i1+1},'5') findstr(NameElectrode{i1+1},'6') findstr(NameElectrode{i1+1},'7') findstr(NameElectrode{i1+1},'8') findstr(NameElectrode{i1+1},'9')])-1);
    if strcmp(tmp1,tmp2)
        PosElectrodeBip(end+1,:)=mean(PosElectrode([i1 i1+1],:));
        NameElectrodeBip{1,end+1}=[NameElectrode{i1+1} NameElectrode{i1}];
        bipole=[bipole [i1+1;i1]];
        NElectrodeBip=NElectrodeBip+1;
    end
end

a_ptsname_pat=fullfile(DirOut,[nameFile '.pts']);
f=fopen(a_ptsname_pat,'w');
fprintf(f,'%s\n','ptsfile');
fprintf(f,'%s\t%s\t%s\n','1','1','1');
fprintf(f,'%s\n',int2str(NElectrode));
for s_c=1:NElectrode
    fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',NameElectrode{s_c},num2str(PosElectrode(s_c,1),2),num2str(PosElectrode(s_c,2),2),num2str(PosElectrode(s_c,3),2),'0','0','0','2','2.0');
end;
fclose(f);

a_ptsname_pat=fullfile(DirOut,[nameFile '_Pos.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrode
    fprintf(f,'%s\t%s\t%s\n',num2str(PosElectrode(s_c,1),2),num2str(PosElectrode(s_c,2),2),num2str(PosElectrode(s_c,3),2));
end;
fclose(f);

a_ptsname_pat=fullfile(DirOut,[nameFile '_Name.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrode
    fprintf(f,'%s\n',NameElectrode{s_c});
end;
fclose(f);

a_ptsname_pat=fullfile(DirOut,[nameFile '_PosBip.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrodeBip
    fprintf(f,'%s\t%s\t%s\n',num2str(PosElectrodeBip(s_c,1),2),num2str(PosElectrodeBip(s_c,2),2),num2str(PosElectrodeBip(s_c,3),2));
end;
fclose(f);

a_ptsname_pat=fullfile(DirOut,[nameFile '_NameBip.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrodeBip
    fprintf(f,'%s\n',NameElectrodeBip{s_c});
end;
fclose(f);

%Compute inverse deformation field
clear matlabbatch
matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def{1}=Deformation;
matlabbatch{1}.spm.util.defs.comp{1}.inv.space{1}=FileSource;
matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname='inverse';
matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr{1}=spm_str_manip(FileSource,'h');
spm_jobman('run',matlabbatch);

%Read deformation field
P=fullfile(spm_str_manip(FileSource,'h'),'y_inverse.nii');
% P=Deformation;
P1=spm_vol([P ',1,1']);
P2=spm_vol([P ',1,2']);
P3=spm_vol([P ',1,3']);
[V1,XYZ]=spm_read_vols(P1);
V2=spm_read_vols(P2);
V3=spm_read_vols(P3);

%Apply tranformation to electrodes
wPosElectrode=PosElectrode;
for i1=1:size(PosElectrode,1)
    D=(XYZ(1,:)-PosElectrode(i1,1)).^2+(XYZ(2,:)-PosElectrode(i1,2)).^2+(XYZ(3,:)-PosElectrode(i1,3)).^2;
    [tmp,order]=sort(D);
    tmp=tmp(1:18);      %cubic neighborhood
    order=order(1:18);
    W=1./tmp;           %weight inverse to distance
    if sum(isinf(W))>0
        W=[1 zeros(1,length(W)-1)];
    end    
    wPosElectrode(i1,:)=[sum(V1(order).*W)./sum(W) sum(V2(order).*W)./sum(W) sum(V3(order).*W)./sum(W)];
end
wPosElectrodeBip=PosElectrodeBip;
for i1=1:size(PosElectrodeBip,1)
    D=(XYZ(1,:)-PosElectrodeBip(i1,1)).^2+(XYZ(2,:)-PosElectrodeBip(i1,2)).^2+(XYZ(3,:)-PosElectrodeBip(i1,3)).^2;
    [tmp,order]=sort(D);
    tmp=tmp(1:18);      %cubic neighborhood
    order=order(1:18);
    W=1./tmp;           %weight inverse to distance
    if sum(isinf(W))>0
        W=[1 zeros(1,length(W)-1)];
    end    
    wPosElectrodeBip(i1,:)=[sum(V1(order).*W)./sum(W) sum(V2(order).*W)./sum(W) sum(V3(order).*W)./sum(W)];
end

a_ptsname_pat=fullfile(DirOut,[nameFile '_MNI.pts']);
f=fopen(a_ptsname_pat,'w');
fprintf(f,'%s\n','ptsfile');
fprintf(f,'%s\t%s\t%s\n','1','1','1');
fprintf(f,'%s\n',int2str(NElectrode));
for s_c=1:NElectrode
    fprintf(f,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',NameElectrode{s_c},num2str(wPosElectrode(s_c,1),2),num2str(wPosElectrode(s_c,2),2),num2str(wPosElectrode(s_c,3),2),'0','0','0','2','2.0');
end;
fclose(f);

a_ptsname_pat=fullfile(DirOut,[nameFile '_Pos_MNI.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrode
    fprintf(f,'%s\t%s\t%s\n',num2str(wPosElectrode(s_c,1),2),num2str(wPosElectrode(s_c,2),2),num2str(wPosElectrode(s_c,3),2));
end;
a_ptsname_pat=fullfile(DirOut,[nameFile '_PosBip_MNI.txt']);
f=fopen(a_ptsname_pat,'w');
for s_c=1:NElectrodeBip
    fprintf(f,'%s\t%s\t%s\n',num2str(wPosElectrodeBip(s_c,1),2),num2str(wPosElectrodeBip(s_c,2),2),num2str(wPosElectrodeBip(s_c,3),2));
end;
fclose(f);

        
