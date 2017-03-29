function ImaGIN_ConnectStim(S)
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
% Copyright (c) 2000-2017 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David

% global defaults
spm('defaults','EEG')

global st

FS1 = spm('FontSize', 14);

% try
%     Prefixe=S.Prefixe;
% catch
%     Prefixe=spm_input('DES feature ', 1,'ERP|Gamma');
%     switch Prefixe
%         case 'ERP'
%             Prefixe='SPM_';
%         case 'Gamma'
%             Prefixe='SPMg_';
%     end
% end
try
    PrefixeLat=S.PrefixeLat;
    PrefixePro=S.PrefixePro;
catch
    PrefixeLat='Latency_';
    PrefixePro='Probability_';
end


try
    PatientList=S.PatientList;
catch
    PatientList=spm_select(inf,'dir');
    S.PatientList=PatientList;
end
Root=fileparts(fileparts(deblank(PatientList(1,:))));

cd(Root)

F = spm_figure('Create','Graphics','Graphics','on');
spm_figure('ColorMap','gray-jet');

draw_anat


WS = spm('WinScale');
uicontrol(F,'Style','Frame','Position',[377 110 180 95].*WS);
uicontrol(F,'Style','Text', 'Position',[382 188 170 016].*WS,'String','Crosshair Position');
uicontrol(F,'Style','Text', 'Position',[382 160 35 025].*WS,'String','mm:');
uicontrol(F,'Style','Text', 'Position',[382 140 35 025].*WS,'String','vx:');
uicontrol(F,'Style','Text', 'Position',[382 118 75 025].*WS,'String','Img Intens.:');
st.mp = uicontrol(F,'Style','edit', 'Position',[418 160 135 025].*WS,'String','','Callback','spm_image(''setposmm'')','ToolTipString','move crosshairs to mm coordinates');
st.vp = uicontrol(F,'Style','edit', 'Position',[418 140 135 025].*WS,'String','','Callback','spm_image(''setposvx'')','ToolTipString','move crosshairs to voxel coordinates');
st.in = uicontrol(F,'Style','Text', 'Position',[428 118  85 025].*WS,'String','');
% st.Callback = 'spm_image(''shopos'');';


h.title=uicontrol(F,'Style','Text', 'Position',[0 720 605 30].*WS,'String','','FontSize',FS1,'BackgroundColor','w');


% Display button
%--------------------------
h.displaybutton=uicontrol(F,'Style','togglebutton','Position',[377 210 82 22].*WS,'String','Build atlas','Callback',@DisplayAtlas,'ToolTipString','Construct DES atlas for current crosshair position');

% Clear button
%--------------------------
h.clearbutton=uicontrol(F,'Style','togglebutton','Position',[463 210 72 22].*WS,'String','Clear atlas','Callback',@ClearAtlas,'ToolTipString','Clear MRI');

% Display N button
%--------------------------
h.Nbutton=uicontrol(F,'Style','togglebutton','Position',[539 210 18 22].*WS,'String','N','Callback',@DisplayNumberStim,'ToolTipString','Display stimulation sites of selected data basis.');

% Latency button
%--------------------------
uicontrol(F,'Style','Frame','Position',[400 400 135 30].*WS);
uicontrol(F,'Style','text','String','Map','HorizontalAlignment','left', ...
    'Position',[404 408 90 19].*WS);
txt_box = cell(6,1);
txt_box{1}='Probability';
txt_box{2}='Latency 0.95';
txt_box{3}='Latency 0.96';
txt_box{4}='Latency 0.97';
txt_box{5}='Latency 0.98';
txt_box{6}='Latency 0.99';
h.latencybutton = uicontrol(F,'Style','popup','String',txt_box, ...
    'Position',[456 408 80 19].*WS, ...
    'Callback',@ChooseLatency,'ToolTipString','Select to display probability maps or median latency maps. Latency estimates depends on the selected proability threshold.');

% Flip button
%--------------------------
uicontrol(F,'Style','Frame','Position',[400 360 135 30].*WS);
uicontrol(F,'Style','text','String','Sum L/R','HorizontalAlignment','left', ...
    'Position',[404 368 90 19].*WS);
txt_box = cell(2,1);
txt_box{1}='No';
txt_box{2}='Yes';
h.flipbutton = uicontrol(F,'Style','popup','String',txt_box, ...
    'Position',[456 368 80 19].*WS, ...
    'Callback',@ChooseFlip,'ToolTipString','Sum both hemispheres to increase number of stimulations around cross-hair.');

% ROI button
%--------------------------
uicontrol(F,'Style','Frame','Position',[400 320 135 30].*WS);
uicontrol(F,'Style','text','String','Radius','HorizontalAlignment','left', ...
    'Position',[404 328 90 19].*WS);
txt_box = cell(6,1);
txt_box{1}='5';
txt_box{2}='10';
txt_box{3}='15';
txt_box{4}='20';
txt_box{5}='25';
txt_box{6}='30';
h.ROIbutton = uicontrol(F,'Style','popup','String',txt_box, ...
    'Position',[456 328 80 19].*WS, 'Value', 3, ...
    'Callback',@ChooseRadius,'ToolTipString','Select radius in mm of ROI sphere around crosshair to count stimulations.');

% Minimal number of stimulations button
%--------------------------
uicontrol(F,'Style','Frame','Position',[400 280 135 30].*WS);
uicontrol(F,'Style','text','String','min(N)','HorizontalAlignment','left', ...
    'Position',[404 288 90 19].*WS);
txt_box = cell(10,1);
txt_box{1}='1';
txt_box{2}='2';
txt_box{3}='3';
txt_box{4}='4';
txt_box{5}='5';
txt_box{6}='6';
txt_box{7}='7';
txt_box{8}='8';
txt_box{9}='9';
txt_box{10}='10';
h.Nminbutton = uicontrol(F,'Style','popup','String',txt_box, ...
    'Position',[456 288 80 19].*WS, 'Value',5,...
    'Callback',@ChooseNmin,'ToolTipString','Select minimal number of recordings for voxels of the atlas.');

% Probability threshold
%--------------------------
uicontrol(F,'Style','Frame','Position',[400 240 135 30].*WS);
uicontrol(F,'Style','text','String','P value','HorizontalAlignment','left', ...
    'Position',[404 248 90 19].*WS);
txt_box = cell(11,1);
txt_box{1}='0';
txt_box{2}='0.1';
txt_box{3}='0.2';
txt_box{4}='0.3';
txt_box{5}='0.4';
txt_box{6}='0.5';
txt_box{7}='0.6';
txt_box{8}='0.7';
txt_box{9}='0.8';
txt_box{10}='0.9';
txt_box{11}='1';
h.pvaluebutton = uicontrol(F,'Style','popup','String',txt_box, ...
    'Position',[456 248 80 19].*WS, 'Value',1,...
    'Callback',@ChoosePvalue,'ToolTipString','Select probability threshold for display.');


% Display button
%--------------------------
h.comparebutton=uicontrol(F,'Style','togglebutton','Position',[377 80 180 19].*WS,'String','Compare to group','Callback',@CompareAtlas,'ToolTipString','Construct patient''s departure from DES atlas');

try
    h.S=S;
end



XYZ=[];
for i1=1:size(PatientList,1)
    directories=dir(deblank(PatientList(i1,:)));
    for i2=3:length(directories)
        tmp=max(findstr(directories(i2).name,'_'));
        if directories(i2).isdir && strcmp(directories(i2).name(1:length(PrefixeLat)),PrefixeLat)
            XYZ(end+1,:)=str2num(directories(i2).name(tmp+1:end));
        end
    end
end


%Initialisation of maps
cd(Root)
n=3;
bb = [[-78 -112 -50];[78 76 85]];
P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
n1=length(bb(1,1):n:bb(2,1));
n2=length(bb(1,2):n:bb(2,2));
n3=length(bb(1,3):n:bb(2,3));
V = fullfile(spm('dir'), 'canonical', 'avg152T1.nii');
V=spm_vol(V);
V.mat=spm_matrix(P);
V.dim=[n1 n2 n3];
V.dt=[16 0];    %float 32
V.fname = 'NumberStim.img';
NImage=zeros(n1,n2,n3);
V=spm_write_vol(V,NImage);

[R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
RCP      = [R(:)';C(:)';P(:)'];
clear R C P
RCP(4,:) = 1;
XYZ1      = V(1).mat(1:3,:)*RCP;
for i1=1:size(XYZ,1)
    D=sqrt(sum((XYZ(i1,:)'*ones(1,size(XYZ1,2))-XYZ1).^2));
    Sel=find(D<=5);
    NImage(Sel)=NImage(Sel)+1;
end
V=spm_vol(V);
V=spm_write_vol(V,NImage);


h.PrefixeLat=PrefixeLat;
h.PrefixePro=PrefixePro;
h.XYZ=XYZ;
h.Root=Root;
h.PatientList=PatientList;
h.st=st;
h.LatencyValues=0.95:0.01:0.99;
guidata(F, h);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_anat
% This function plots anatomical MRI

global st

spm_orthviews('Reset');
st.callback = 'spm_image(''shopos'');';
P = fullfile(spm('dir'),'canonical','avg152T1.nii');
spms = spm_vol(P);
spm_orthviews('Image',spms);
%     spm_orthviews('Image',spms,[0.05 0.05 0.9 0.45]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChoosePvalue(hObject, events)

global st

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

W=h.W;
DataBin2=h.DataBin2;
xSPM=h.xSPM;

Latency=h.Latency;

z=get(h.pvaluebutton);
Zmax=str2num(z.String{z.Value});
h.Zmax=Zmax;    
    
%Display image
    i=find(DataBin2>=Zmax);
    xSPM.Z=W(i)';
    [i,j,k]=ind2sub(size(W),i);
    xSPM.XYZ=[i j k]';
    xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];    
    



CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;


h.xSPM=xSPM;
guidata(F, h);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChooseNmin(hObject, events)

global st

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

z=get(h.Nminbutton);
Nmin=str2num(z.String{z.Value});
h.Nmin=Nmin;

SelectionMRI=h.SelectionMRI;
Symmetry=h.Symmetry;
Latency=h.Latency;
Dmax=h.Dmax;
Zmax=h.Zmax;
DataBin=h.DataBin;
DataLat=h.DataLat;
N=h.N;
NLat=h.NLat;
xSPM=h.xSPM;
LatencyValues=h.LatencyValues;
    
    N=h.N;
    Dmax=h.Dmax;
    Symmetry=h.Symmetry;
    DataBin2=h.DataBin;
    Root=h.Root;
    SelectionMRI=h.SelectionMRI;
    n1=h.n1;
    n2=h.n2;
    n3=h.n3;
    
    
%Cut off in voxels with too few recordings
DataBin2=DataBin;
DataBin2(find(N<Nmin))=NaN;
%Show target
DataBin2(SelectionMRI)=1;
DataLat2=DataLat;
for i1=1:length(DataLat2)
    DataLat2{i1}(find(NLat{i1}<Nmin))=NaN;
    DataLat2{i1}(SelectionMRI)=0;
end


%save images
if ~Symmetry
    V.fname = ['Probability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['Probability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V.fname);
spm_write_vol(V,DataBin2);
%smooth image
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'EEG');
spm_jobman('run', matlabbatch);
if ~Symmetry
    V.fname = ['sProbability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['sProbability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V.fname);
DataBin2=spm_read_vols(V);

for i1=1:length(LatencyValues)
    if ~Symmetry
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V);
    V=spm_write_vol(V,DataLat2{i1});
    %smooth image
    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm('defaults', 'EEG');
    spm_jobman('run', matlabbatch);
    if ~Symmetry
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V.fname);
    DataLat2{i1}=spm_read_vols(V);
end

%Display image
switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        I2=DataLat2{tmp};
        Zmax=-1;
    case 'Pro'
        I2=DataBin2;
end
W=I2;
xSPM.M          = V.mat;
    
    i=find(~isnan(W));
    xSPM.Z=W(i)';
    [i,j,k]=ind2sub(size(W),i);
    xSPM.XYZ=[i j k]';
    xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];    
    



CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;


h.xSPM=xSPM;
h.W=W;
h.DataBin2=DataBin2;
h.DataLat2=DataLat2;
h.N=N;
h.NLat=NLat;
guidata(F, h);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChooseLatency(hObject, events)

global st

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

PrefixeLat=h.PrefixeLat;
PrefixePro=h.PrefixePro;
z=get(h.latencybutton);
Latency=z.String{z.Value};
switch Latency
    case 'Latency 0.95'
        Prefixe=[PrefixeLat '950_'];
    case 'Latency 0.96'
        Prefixe=[PrefixeLat '960_'];
    case 'Latency 0.97'
        Prefixe=[PrefixeLat '970_'];
    case 'Latency 0.98'
        Prefixe=[PrefixeLat '980_'];
    case 'Latency 0.99'
        Prefixe=[PrefixeLat '990_'];
    case 'Probability'
        Prefixe=PrefixePro;
end
h.Latency=Latency;

Root=h.Root;
PatientList=h.PatientList;
DataBase=h.DataBase;
Selection=h.Selection;
Target=h.Target;
Nmin=h.Nmin;
Symmetry=h.Symmetry;
SelectionMRI=h.SelectionMRI;
Zmax=h.Zmax;
Dmax=h.Dmax;
xSPM=h.xSPM;
n1=h.n1;
n2=h.n2;
n3=h.n3;
LatencyValues=h.LatencyValues;



    
    
    %Display image
switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        I2=h.DataLat2{tmp};
        Zmax=-1;
    case 'Pro'
        I2=h.DataBin2;
end
    W=I2;
    i=find(W>Zmax);
    xSPM.Z=W(i)';
    [i,j,k]=ind2sub(size(W),i);
    xSPM.XYZ=[i j k]';
    xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];

    



CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;

switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        set(h.title,'String',['Median latency (p=' num2str(LatencyValues(tmp)) ') at ' sprintf('%d %d %d',round(Target))]);
    case 'Pro'
        set(h.title,'String',['Probability at ' sprintf('%d %d %d',round(Target))]);
end
set(h.displaybutton,'Value',0);


h.xSPM=xSPM;
h.W=W;
guidata(F, h);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChooseFlip(hObject, events)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

z=get(h.flipbutton);
Symmetry=z.String{z.Value};
switch Symmetry
    case 'Yes'
        Symmetry=1;
    case 'No'
        Symmetry=0;
end
h.Symmetry=Symmetry;
guidata(F, h);
ChooseRadius


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChooseRadius(hObject, events)

global st

F  = spm_figure('GetWin','Graphics');
h = guidata(F);

z=get(h.ROIbutton);
Dmax=str2num(z.String{z.Value});
h.Dmax=Dmax;

Root=h.Root;
Target=h.Target;
DataBase=h.DataBase;
Nmin=h.Nmin;
Symmetry=h.Symmetry;
SelectionMRI=h.SelectionMRI;
Latency=h.Latency;
Zmax=h.Zmax;
xSPM=h.xSPM;
n1=h.n1;
n2=h.n2;
n3=h.n3;
ProbabilityThreshold=h.ProbabilityThreshold;
LatencyValues=h.LatencyValues;


D=zeros(1,size(DataBase,2));
for i1=1:size(DataBase,2)
    tmp=DataBase(1,i1).xyz;
    if Symmetry && sign(Target(1))~=sign(DataBase(1,i1).xyz(1))
        tmp(1)=-tmp(1);
    end
    D(i1)=sqrt(sum((Target-tmp).^2));
end
Selection=find(D<=Dmax);


%Initialisation of maps
cd(Root)
n=3;
bb = [[-78 -112 -50];[78 76 85]];
P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
n1=length(bb(1,1):n:bb(2,1));
n2=length(bb(1,2):n:bb(2,2));
n3=length(bb(1,3):n:bb(2,3));
h.n1=n1;
h.n2=n2;
h.n3=n3;


V = fullfile(spm('dir'), 'canonical', 'avg152T1.nii');
V=spm_vol(V);
V.mat=spm_matrix(P);
V.dim=[n1 n2 n3];
V.dt=[16 0];    %float 32
for i1=1:length(LatencyValues)
    if ~Symmetry
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    spm_write_vol(V,zeros(n1,n2,n3));
end
if ~Symmetry
    V.fname = ['Probability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['Probability_Distance' num2str(Dmax) '_Sym.img'];
end
spm_write_vol(V,zeros(n1,n2,n3));


N=zeros(n1,n2,n3);
NLat=zeros(n1,n2,n3);
DataBin=zeros(n1,n2,n3);
clear DataLat NLat
for i1=1:length(LatencyValues)
    DataLat{i1}=zeros(n1,n2,n3);
    NLat{i1}=zeros(n1,n2,n3);
end
clear tmp
for i0=1:length(Selection)
    if sign(Target(1))==sign(DataBase(Selection(i0)).xyz(1))
        P1=spm_vol(fullfile(DataBase(Selection(i0)).fnamepro,'sample_0.img'));
    else
        P1=spm_vol(fullfile(DataBase(Selection(i0)).fnamepro,'fsample_0.img'));
    end
    [Vb1,XYZ1] = spm_read_vols(P1);
    if i0==1
        D=sqrt(sum((Target'*ones(1,size(XYZ1,2))-XYZ1).^2));
        %             SelectionMRI=find(D<=Dmax);
        SelectionMRI=find(D<=5);
        h.SelectionMRI=SelectionMRI;
    end
    % Activated voxels
    Q1=find(Vb1>=ProbabilityThreshold);
    % Recorded voxels
    Q3=find(~isnan(Vb1));
    tmpBin=zeros(n1,n2,n3);
    tmpBin(Q1)=1;
    tmpNBin=zeros(n1,n2,n3);
    tmpNBin(Q3)=1;
    N=N+tmpNBin;
    DataBin=DataBin+tmpBin;
    for i1=1:length(LatencyValues)
        if sign(Target(1))==sign(DataBase(Selection(i0)).xyz(1))
            P2=spm_vol(fullfile(DataBase(Selection(i0)).fnamelat{i1},'ssample_0.img'));
        else
            P2=spm_vol(fullfile(DataBase(Selection(i0)).fnamelat{i1},'fssample_0.img'));
        end
        Vb2 = spm_read_vols(P2);
        QLat=find(~isnan(Vb2));
        tmpNLatBin=zeros(n1,n2,n3);
        tmpNLatBin(QLat)=1;
        tmpLat=zeros(n1,n2,n3);
        tmpLat(QLat)=Vb2(QLat);
        tmpLat(isnan(tmpLat))=0;
        DataLat{i1}=DataLat{i1}+tmpLat;
        NLat{i1}=NLat{i1}+tmpNLatBin;
    end
    
        
    
    
end


%Normalisation
DataBin=DataBin./N;
%Cut off in voxels with too few recordings
DataBin2=DataBin;
DataBin2(find(N<Nmin))=NaN;
%Show target
DataBin2(SelectionMRI)=1;
for i1=1:length(LatencyValues)
    DataLat{i1}=DataLat{i1}./NLat{i1};
    DataLat2{i1}=DataLat{i1};
    DataLat2{i1}(find(NLat{i1}<Nmin))=NaN;
    DataLat2{i1}(SelectionMRI)=0;
end

%save images
if ~Symmetry
    V.fname = ['Probability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['Probability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V);
spm_write_vol(V,DataBin2);
%smooth image
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'EEG');
spm_jobman('run', matlabbatch);
if ~Symmetry
    V.fname = ['sProbability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['sProbability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V.fname);
DataBin2=spm_read_vols(V);

for i1=1:length(LatencyValues)
    if ~Symmetry
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V);
    V=spm_write_vol(V,DataLat2{i1});
    %smooth image
    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm('defaults', 'EEG');
    spm_jobman('run', matlabbatch);
    if ~Symmetry
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V.fname);
    DataLat2{i1}=spm_read_vols(V);
end

%Display image
switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        I2=DataLat2{tmp};
        Zmax=-1;
    case 'Pro'
        I2=DataBin2;
end
W=I2;
xSPM.M          = V.mat;
    
    xSPM.Z=zeros(1,0);
    xSPM.XYZ=zeros(3,0);
    xSPM.XYZmm=zeros(3,0);
    i=find(W>Zmax);
    xSPM.Z=W(i)';
    [i,j,k]=ind2sub(size(W),i);
    xSPM.XYZ=[i j k]';
    xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    





CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;

switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        set(h.title,'String',['Median latency (p=' num2str(LatencyValues(tmp)) ') at ' sprintf('%d %d %d',round(Target))]);
    case 'Pro'
        set(h.title,'String',['Probability at ' sprintf('%d %d %d',round(Target))]);
end
set(h.displaybutton,'Value',0);


h.xSPM=xSPM;
h.W=W;
h.DataBin=DataBin;
h.DataLat=DataLat;
h.N=N;
h.NLat=NLat;
guidata(F, h);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ClearAtlas(hObject, events)

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
draw_anat;
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;
set(h.title,'String','');
set(h.clearbutton,'Value',0);
guidata(F, h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisplayNumberStim(hObject, events)

global st

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
Root=h.Root;
cd(Root)



V = 'NumberStim.img';
M=spm_vol(V);
I=spm_read_vols(M);

W=I;
xSPM.M          = M.mat;

i=find(W>0);
xSPM.Z=W(i)';
[i,j,k]=ind2sub(size(W),i);
xSPM.XYZ=[i j k]';
xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];



draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;

set(h.title,'String','');
set(h.Nbutton,'Value',0);

h.xSPM=xSPM;
h.W=W;
guidata(F, h);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DisplayAtlas(hObject, events)

global st defaults

Target=st.centre';

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
PrefixeLat=h.PrefixeLat;
PrefixePro=h.PrefixePro;
h.Target=Target;
% st=h.st;
PatientList=h.PatientList;
Root=h.Root;
LatencyValues=h.LatencyValues;
ProbabilityThreshold=0.999;

z=get(h.latencybutton);
Latency=z.String{z.Value};
switch Latency
    case 'Latency 0.95'
        Prefixe=[PrefixeLat '950_'];
    case 'Latency 0.96'
        Prefixe=[PrefixeLat '960_'];
    case 'Latency 0.97'
        Prefixe=[PrefixeLat '970_'];
    case 'Latency 0.98'
        Prefixe=[PrefixeLat '980_'];
    case 'Latency 0.99'
        Prefixe=[PrefixeLat '990_'];
    case 'Probability'
        Prefixe=PrefixePro;
end
h.Latency=Latency;

z=get(h.flipbutton);
Symmetry=z.String{z.Value};
switch Symmetry
    case 'Yes'
        Symmetry=1;
    case 'No'
        Symmetry=0;
end
h.Symmetry=Symmetry;

z=get(h.ROIbutton);
Dmax=str2num(z.String{z.Value});
h.Dmax=Dmax;

z=get(h.Nminbutton);
Nmin=str2num(z.String{z.Value});
h.Nmin=Nmin;

z=get(h.pvaluebutton);
Zmax=str2num(z.String{z.Value});
h.Zmax=Zmax;


clear DataBase
n=0;
for i1=1:size(PatientList,1)
    directories=dir(deblank(PatientList(i1,:)));
    for i2=3:length(directories)
        tmp=max(findstr(directories(i2).name,'_'));
        if directories(i2).isdir && strcmp(directories(i2).name(1:length(Prefixe)),[Prefixe])
            n=n+1;
            for i3=1:length(LatencyValues)
                DataBase(1,n).fnamelat{i3}=fullfile(deblank(PatientList(i1,:)),[PrefixeLat num2str(1000*LatencyValues(i3)) '_' directories(i2).name(length(Prefixe)+1:end)]);
            end
            DataBase(1,n).fnamepro=fullfile(deblank(PatientList(i1,:)),[PrefixePro directories(i2).name(length(Prefixe)+1:end)]);
            DataBase(1,n).xyz=str2num(directories(i2).name(tmp+1:end));
        end
    end
end
h.DataBase=DataBase;

D=zeros(1,size(DataBase,2));
for i1=1:size(DataBase,2)
    tmp=DataBase(1,i1).xyz;
    if Symmetry && sign(Target(1))~=sign(DataBase(1,i1).xyz(1))
        tmp(1)=-tmp(1);
    end
    D(i1)=sqrt(sum((Target-tmp).^2));
end
Selection=find(D<=Dmax);
h.Selection=Selection;


%Initialisation of maps
cd(Root)
n=3;
bb = [[-78 -112 -50];[78 76 85]];
P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
n1=length(bb(1,1):n:bb(2,1));
n2=length(bb(1,2):n:bb(2,2));
n3=length(bb(1,3):n:bb(2,3));
h.n1=n1;
h.n2=n2;
h.n3=n3;


V = fullfile(spm('dir'), 'canonical', 'avg152T1.nii');
V=spm_vol(V);
V.mat=spm_matrix(P);
V.dim=[n1 n2 n3];
V.dt=[16 0];    %float 32
for i1=1:length(LatencyValues)
    if ~Symmetry
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    spm_write_vol(V,zeros(n1,n2,n3));
end
if ~Symmetry
    V.fname = ['Probability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['Probability_Distance' num2str(Dmax) '_Sym.img'];
end
spm_write_vol(V,zeros(n1,n2,n3));


N=zeros(n1,n2,n3);
NLat=zeros(n1,n2,n3);
DataBin=zeros(n1,n2,n3);
clear DataLat NLat
for i1=1:length(LatencyValues)
    DataLat{i1}=zeros(n1,n2,n3);
    NLat{i1}=zeros(n1,n2,n3);
end
clear tmp
for i0=1:length(Selection)
    if sign(Target(1))==sign(DataBase(Selection(i0)).xyz(1))
        P1=spm_vol(fullfile(DataBase(Selection(i0)).fnamepro,'sample_0.img'));
    else
        P1=spm_vol(fullfile(DataBase(Selection(i0)).fnamepro,'fsample_0.img'));
    end
    [Vb1,XYZ1] = spm_read_vols(P1);
    if i0==1
        D=sqrt(sum((Target'*ones(1,size(XYZ1,2))-XYZ1).^2));
        %             SelectionMRI=find(D<=Dmax);
        SelectionMRI=find(D<=5);
        h.SelectionMRI=SelectionMRI;
    end
    % Activated voxels
    Q1=find(Vb1>=ProbabilityThreshold);
    % Recorded voxels
    Q3=find(~isnan(Vb1));
    tmpBin=zeros(n1,n2,n3);
    tmpBin(Q1)=1;
    tmpNBin=zeros(n1,n2,n3);
    tmpNBin(Q3)=1;
    N=N+tmpNBin;
    DataBin=DataBin+tmpBin;
    for i1=1:length(LatencyValues)
        if sign(Target(1))==sign(DataBase(Selection(i0)).xyz(1))
            P2=spm_vol(fullfile(DataBase(Selection(i0)).fnamelat{i1},'ssample_0.img'));
        else
            P2=spm_vol(fullfile(DataBase(Selection(i0)).fnamelat{i1},'fssample_0.img'));
        end
        Vb2 = spm_read_vols(P2);
        QLat=find(~isnan(Vb2));
        tmpNLatBin=zeros(n1,n2,n3);
        tmpNLatBin(QLat)=1;
        tmpLat=zeros(n1,n2,n3);
        tmpLat(QLat)=Vb2(QLat);
        tmpLat(isnan(tmpLat))=0;
        DataLat{i1}=DataLat{i1}+tmpLat;
        NLat{i1}=NLat{i1}+tmpNLatBin;
    end
    
        
    
    
end


%Normalisation
DataBin=DataBin./N;
%Cut off in voxels with too few recordings
DataBin2=DataBin;
DataBin2(find(N<Nmin))=NaN;
%Show target
DataBin2(SelectionMRI)=1;
for i1=1:length(LatencyValues)
    DataLat{i1}=DataLat{i1}./NLat{i1};
    DataLat2{i1}=DataLat{i1};
    DataLat2{i1}(find(NLat{i1}<Nmin))=NaN;
    DataLat2{i1}(SelectionMRI)=0;
end

%save images
if ~Symmetry
    V.fname = ['Probability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['Probability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V);
spm_write_vol(V,DataBin2);
%smooth image
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'EEG');
spm_jobman('run', matlabbatch);
if ~Symmetry
    V.fname = ['sProbability_Distance' num2str(Dmax) '.img'];
else
    V.fname = ['sProbability_Distance' num2str(Dmax) '_Sym.img'];
end
V=spm_vol(V.fname);
DataBin2=spm_read_vols(V);

for i1=1:length(LatencyValues)
    if ~Symmetry
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['Latency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V);
    V=spm_write_vol(V,DataLat2{i1});
    %smooth image
    clear matlabbatch
    matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 1;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm('defaults', 'EEG');
    spm_jobman('run', matlabbatch);
    if ~Symmetry
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '.img'];
    else
        V.fname = ['sLatency_' num2str(1000*LatencyValues(i1)) '_Distance' num2str(Dmax) '_Sym.img'];
    end
    V=spm_vol(V.fname);
    DataLat2{i1}=spm_read_vols(V);
end

%Display image
switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        I2=DataLat2{tmp};
        Zmax=-1;
    case 'Pro'
        I2=DataBin2;
end
W=I2;
xSPM.M          = V.mat;

xSPM.Z=zeros(1,0);
xSPM.XYZ=zeros(3,0);
xSPM.XYZmm=zeros(3,0);
i=find(W>Zmax);
xSPM.Z=W(i)';
[i,j,k]=ind2sub(size(W),i);
xSPM.XYZ=[i j k]';
xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];

CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;

switch Latency(1:3)
    case 'Lat'
        tmp=find(LatencyValues==str2double(Latency(end-3:end)));
        set(h.title,'String',['Median latency (p=' num2str(LatencyValues(tmp)) ') at ' sprintf('%d %d %d',round(Target))]);
    case 'Pro'
        set(h.title,'String',['Probability at ' sprintf('%d %d %d',round(Target))]);
end
set(h.displaybutton,'Value',0);

h.xSPM=xSPM;
h.W=W;
h.DataBin=DataBin;
h.DataLat=DataLat;
h.DataBin2=DataBin2;
h.DataLat2=DataLat2;
h.N=N;
h.NLat=NLat;
h.ProbabilityThreshold=ProbabilityThreshold;
guidata(F, h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CompareAtlas(hObject, events)

global st defaults

F  = spm_figure('GetWin','Graphics');
h = guidata(F);
Prefixe=h.PrefixePro;
% st=h.st;
PatientList=h.PatientList;
DataBase=h.DataBase;
Root=h.Root;
LatencyValues=h.LatencyValues;
ProbabilityThreshold=h.ProbabilityThreshold;

PrefixeLat=h.PrefixeLat;
PrefixePro=h.PrefixePro;
z=get(h.latencybutton);
Latency=z.String{z.Value};
% switch Latency(1:3)
%     case 'Lat'
%         tmp=find(LatencyValues==str2double(Latency(end-3:end)));
%         Prefixe=[PrefixeLat num2str(1000*LatencyValues(tmp)) '_'];
%     case 'Pro'
%         Prefixe=PrefixePro;
% end
% 
% 
% Prefixe=PrefixePro;


h.Latency=Latency;


z=get(h.flipbutton);
Symmetry=z.String{z.Value};
switch Symmetry
    case 'Yes'
        Symmetry=1;
    case 'No'
        Symmetry=0;
end
h.Symmetry=Symmetry;

z=get(h.ROIbutton);
Dmax=str2num(z.String{z.Value});
h.Dmax=Dmax;

z=get(h.Nminbutton);
Nmin=str2num(z.String{z.Value});
h.Nmin=Nmin;

z=get(h.pvaluebutton);
Zmax=str2num(z.String{z.Value});
h.Zmax=Zmax;




%Select patient of interest
Patient=spm_select(1,'dir');
PatientName=spm_str_manip(spm_str_manip(Patient,'h'),'t');

clear DataBase2
n=0;
for i1=1:size(Patient,1)
    directories=dir(deblank(Patient(i1,:)));
    for i2=3:length(directories)
        tmp=max(findstr(directories(i2).name,'_'));
        if directories(i2).isdir && strcmp(directories(i2).name(1:length(Prefixe)),[Prefixe])
            n=n+1;
            DataBase2(1,n).fname=fullfile(deblank(Patient(i1,:)),[Prefixe directories(i2).name(length(Prefixe)+1:end)]);
            DataBase2(1,n).xyz=str2num(directories(i2).name(tmp+1:end));
        end
    end
end
h.DataBase2=DataBase2;

Select=[];
for i1=1:size(PatientList,1)
    if strcmp(deblank(PatientList(i1,:)),deblank(Patient))
        Select=i1;
    end
end
Select=setdiff(1:size(PatientList,1),Select);

clear DataBase3
n=0;
for i1=1:size(PatientList(Select,:),1)
    directories=dir(deblank(PatientList(Select(i1),:)));
    for i2=3:length(directories)
        tmp=max(findstr(directories(i2).name,'_'));
        if directories(i2).isdir && strcmp(directories(i2).name(1:length(Prefixe)),[Prefixe])
            n=n+1;
            DataBase3(1,n).fname=fullfile(deblank(PatientList(Select(i1),:)),[Prefixe directories(i2).name(length(Prefixe)+1:end)]);
            DataBase3(1,n).xyz=str2num(directories(i2).name(tmp+1:end));
        end
    end
end
h.DataBase3=DataBase3;


%Initialisation of patient's atlas map
cd(Root)
n=3;
bb = [[-78 -112 -50];[78 76 85]];
P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
n1=length(bb(1,1):n:bb(2,1));
n2=length(bb(1,2):n:bb(2,2));
n3=length(bb(1,3):n:bb(2,3));
h.n1=n1;
h.n2=n2;
h.n3=n3;


V = fullfile(spm('dir'), 'canonical', 'avg152T1.nii');
V=spm_vol(V);
V.mat=spm_matrix(P);
V.dim=[n1 n2 n3];
V.dt=[16 0];    %float 32
if ~Symmetry
    V.fname = ['Similarity_' PatientName '.img'];
else
    V.fname = ['Similarity_' PatientName '_Sym.img'];
end
V=spm_write_vol(V,zeros(n1,n2,n3));

Similarity=zeros(length(DataBase2),1);
ElectrodePos=zeros(length(DataBase2),3);
for i00=1:size(DataBase2,2)
    
    
    Target=DataBase2(1,i00).xyz;
    ElectrodePos(i00,:)=Target;
    
    
    D=zeros(1,size(DataBase3,2));
    for i1=1:size(DataBase3,2)
        tmp=DataBase3(1,i1).xyz;
        if Symmetry && sign(Target(1))~=sign(DataBase3(1,i1).xyz(1))
            tmp(1)=-tmp(1);
        end
        D(i1)=sqrt(sum((Target-tmp).^2));
    end
    Selection=find(D<=Dmax);
    
    
    for i000=1:size(DataBase2,1)
        
        disp([i000 size(DataBase2,1) i00 size(DataBase2,2)])
        
        %%%%%%%%%%%%%%%%%%%%%
        %Read patient's connectivity
        %%%%%%%%%%%%%%%%%%%%%
        
        
        if sign(Target(1))==sign(DataBase2(i000,i00).xyz(1))
            P1=spm_vol(fullfile(DataBase2(i000,i00).fname,'sample_0.img'));
        else
            P1=spm_vol(fullfile(DataBase2(i000,i00).fname,'fsample_0.img'));
        end
        [Vb1,XYZ1] = spm_read_vols(P1);
        % Activated voxels
        Q10=find(Vb1>=ProbabilityThreshold);
        % Recorded voxels
        Q30=find(~isnan(Vb1));
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%
        %Read atlas connectivity
        %%%%%%%%%%%%%%%%%%%%%
        
        
        N=zeros(n1,n2,n3);
        DataBin=zeros(n1,n2,n3);
        for i0=1:length(Selection)
            if sign(Target(1))==sign(DataBase3(i000,Selection(i0)).xyz(1))
                P1=spm_vol(fullfile(DataBase3(i000,Selection(i0)).fname,'sample_0.img'));
            else
                P1=spm_vol(fullfile(DataBase3(i000,Selection(i0)).fname,'fsample_0.img'));
            end
            [Vb1,XYZ1] = spm_read_vols(P1);
            
            % Activated voxels
            Q1=find(Vb1>=ProbabilityThreshold);
            % Recorded voxels
            Q3=find(~isnan(Vb1));
            
            tmpBin=zeros(n1,n2,n3);
            tmpBin(Q1)=1;
            tmpNBin=zeros(n1,n2,n3);
            tmpNBin(Q3)=1;
            N=N+tmpNBin;
            DataBin=DataBin+tmpBin;
            
        end
        
        %Normalisation
        DataBin=DataBin./N;
        
        %Cut off in voxels with too few recordings
        DataBin2=DataBin;
        DataBin2(find(N<Nmin))=NaN;
        
        %Read activated voxels of the atlas
        Q1=find(DataBin2>=Zmax);
        
        
        %inclusive mask with patient's voxel
        Q11=intersect(Q1,Q30);
        
        %Measure of similarity
        Simi=length(intersect(Q11,Q10))/length(Q11);
        
        Similarity(i00,i000)=1-Simi;
        
    end
end
Similarity=mean(Similarity,2);

%write mat file
D = [];
D.Atlas='Human';
D.descrip.date = ' / / ';
D.descrip.time = ' / / ';
Nchannels = length(Similarity);
D.Fsample=1;
Nsamples = 1; 
D.channels = repmat(struct('bad', 0), 1, 1);
for i2=1:length(Similarity)
    D.channels(i2).label=num2str(i2);
end
D.channelOrder = 1:Nchannels;
D.data.fnamedat = ['Similarity_' num2str(10*Zmax) '_' PatientName '.dat'];
D.data.datatype = 'float32-le';
D.Nchannels = Nchannels;
D.Nsamples = Nsamples;
nsampl=Nsamples;
nchan=D.Nchannels;
D.trials.label = 'Undefined';
D.trials.events = [];
D.trials.onset = 1/D.Fsample;
%Continuous
datafile = file_array(D.data.fnamedat, [nchan nsampl], 'float32-le');
% physically initialise file
datafile(end,end) = 0;
offset = 1;
nblocksamples = size(Similarity,2);
datafile(:, offset:(offset+nblocksamples-1)) = Similarity;
offset = offset+nblocksamples;


%Electrodes
D.sensors.eeg.pnt=ElectrodePos;
for i1=1:length(D.channels)
    D.sensors.eeg.label{1,i1}=D.channels(i1).label;
end
D.sensors.eeg.unit='mm';

D.Atlas='Human';

%--------- Create meeg object
D.fname = ['Similarity_' num2str(10*Zmax) '_' PatientName '.mat'];

D = meeg(D);
save(D)    

clear SS
SS.n=3;
SS.TimeWindow=0;
SS.TimeWindowWidth=0;
SS.interpolate_bad=1;
SS.SizeSphere=5;
SS.SizeHorizon=10;
SS.Atlas='Human';
if isdir(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName]))
    cd(Root)
    rmdir(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName]),'s')
end
SS.Fname=fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName]);
ImaGIN_spm_eeg_convertmat2ana_3D(SS)

if ~Symmetry
    copyfile(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName],'sample_0.hdr'),fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName '.hdr']))
    copyfile(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName],'sample_0.img'),fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName '.img']))
else
    copyfile(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName],'sample_0.hdr'),fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName '_Sym.hdr']))
    copyfile(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName],'sample_0.img'),fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName '_Sym.img']))
end
rmdir(fullfile(Root,['Similarity_' num2str(10*Zmax) '_' PatientName]),'s')


if ~Symmetry
    V.fname = ['Similarity_' num2str(10*Zmax) '_' PatientName '.img'];
else
    V.fname = ['Similarity_' num2str(10*Zmax) '_' PatientName '_Sym.img'];
end
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm('defaults', 'EEG');
spm_jobman('run', matlabbatch);
if ~Symmetry
    V.fname = ['sSimilarity_' num2str(10*Zmax) '_' PatientName '.img'];
else
    V.fname = ['sSimilarity_' num2str(10*Zmax) '_' PatientName '_Sym.img'];
end
V=spm_vol(V.fname);
I2=spm_read_vols(V);

%Display image
W=I2;
xSPM.M          = V.mat;

xSPM.Z=zeros(1,0);
xSPM.XYZ=zeros(3,0);
xSPM.XYZmm=zeros(3,0);
i=find(W>-inf);
xSPM.Z=W(i)';
[i,j,k]=ind2sub(size(W),i);
xSPM.XYZ=[i j k]';
xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];


CurrentTarget=st.centre;
spm_orthviews('Reset');
draw_anat;
spm_orthviews('addblobs',1,xSPM.XYZ,xSPM.Z,xSPM.M);
spm_orthviews('Redraw');
spm_orthviews('Reposition',CurrentTarget)
st.mp=h.st.mp;
st.vp=h.st.vp;
st.in=h.st.in;


set(h.title,'String',['Similarity ' num2str(10*Zmax) '_' PatientName]);
set(h.comparebutton,'Value',0);


h.xSPM=xSPM;
h.W=W;
h.DataBin=DataBin;
h.N=N;
guidata(F, h);


