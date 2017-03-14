function ImaGIN_Disp3D

Fgraph  = spm_figure('GetWin','Graphics');
Finter  = spm_figure('GetWin','Interactive');
spm_figure('Clear',Finter)

tmp=spm_input('Image of p-values ',1,'Yes|No');
if strcmp(tmp,'Yes')
    FlagStat=1;
else
    FlagStat=0;
end

V = spm_select(1, 'image', 'Select 3D data file');

tmp=spm_str_manip(V,'h');
cd(tmp)

tmp=spm_input('Select another 3D data file ',1,'Yes|No');
if strcmp(tmp,'Yes')
    V2 = spm_select(1, 'image', 'Select 3D data file');
end

%check if synchrony
tmp=spm_input('Connections associated ',1,'Yes|No');
if strcmp(tmp,'Yes')
    FlagSyn=1;
    FileSyn = spm_select(1, '.syn', 'Select synchrony file (.syn)');
else
    FlagSyn=0;
end
if FlagSyn
    Syn=load(FileSyn,'-mat');
    if FlagStat
        u  = spm_input(['threshold sy {',Syn.STAT,' or p value}'],'+0','r',0.001,1);
        %Bonferonni
        u=u/length(Syn.S1);
        if u <= 1; u = spm_u(u^(1/Syn.n),Syn.df,Syn.STAT); end
        PairPos=find(Syn.S1>u);
        PairNeg=find(Syn.S2>u);
    else
        mas=max(Syn.S1(:));
        mis=max(Syn.S2(:));
        Zmaxs = str2num(spm_input(sprintf('+ threshold sy (max=%0.3f) ',mas), '+1', 's',0));
        Zmins = str2num(spm_input(sprintf('- threshold sy (min=%0.3f) ',mis), '+1', 's',0));
        PairPos=find(Syn.S1>Zmaxs);
        PairNeg=find(Syn.S2>Zmins);
    end
    M=ImaGIN_ConnectivityMatrix(size(Syn.Pos,1));
end

%plot electrodes
tmp=spm_input('Plot electrodes ',1,'Yes|No');
if strcmp(tmp,'Yes')
    FlagElec=1;
    FileElec = spm_select(inf, '.txt', 'Select electrode positions file (.txt)');
else
    FlagElec=0;
end



V=spm_vol(V);
xSPM.M          = V.mat;
DIM        = V.dim';
W=spm_read_vols(V(1));
% [i,j,k]=find(isfinite(W));
if exist('V2')
    V2=spm_vol(V2);
    W2=spm_read_vols(V2(1));
end
if FlagStat
    if FlagSyn
        u  = spm_input(['threshold im {',Syn.STAT,' or p value}'],'+0','r',0.001,1);
        if u <= 1; u = spm_u(u^(1/Syn.n),Syn.df,Syn.STAT); end
    else
        STAT='T';
        n=1;
        clear df
        df(1)=1;
        df(2)=spm_input('Degrees of freedom','+1', 'r');
        u  = spm_input(['threshold im {','T',' or p value}'],'+0','r',0.001,1);
        if u <= 1; u = spm_u(u^(1/n),df,STAT); end
    end
    i1=find(W>u);
    if exist('W2')
        i2=find(W2>u);
        i=[i1;i2];
        xSPM1=xSPM;
        xSPM2=xSPM;
        xSPM.Z=[W(i1)' W2(i2)'];
        [i,j,k]=ind2sub(V.dim,i);
        xSPM.XYZ=[i j k]';
        xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
        xSPM1.Z=W(i1)';
        [i,j,k]=ind2sub(V.dim,i1);
        xSPM1.XYZ=[i j k]';
        xSPM1.XYZmm=xSPM1.M(1:3,:)*[xSPM1.XYZ; ones(1,size(xSPM1.XYZ,2))];
        xSPM2.Z=W2(i2)';
        [i,j,k]=ind2sub(V.dim,i2);
        xSPM2.XYZ=[i j k]';
        xSPM2.XYZmm=xSPM2.M(1:3,:)*[xSPM2.XYZ; ones(1,size(xSPM2.XYZ,2))];
    else
        i=i1;
        xSPM.Z=W(i)';
        [i,j,k]=ind2sub(V.dim,i);
        xSPM.XYZ=[i j k]';
        xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    end
else
    ma=max(W(:));
    mi=min(W(:));
    Zmax=0;
    Zmin=0;
    if ma>0
        Zmax = str2num(spm_input(sprintf('+ threshold im (max=%0.3f) ',ma), '+1', 's'));
    end
    if mi<0
        Zmin = str2num(spm_input(sprintf('- threshold im (min=%0.3f) ',mi), '+1', 's'));
    end
    xSPM.Z=zeros(1,0);
    xSPM.XYZ=zeros(3,0);
    xSPM.XYZmm=zeros(3,0);
    if ma>Zmax&mi>=Zmin
        i=find(W>Zmax);
        xSPM.Z=W(i)';
        [i,j,k]=ind2sub(V.dim,i);
        xSPM.XYZ=[i j k]';
        xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    elseif ma<=Zmax&mi<Zmin
        i=find(W<Zmin);
        xSPM.Z=W(i)';
        [i,j,k]=ind2sub(V.dim,i);
        xSPM.XYZ=[i j k]';
        xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    elseif ma>Zmax&mi<Zmin
        i1=find(W>Zmax);
        i2=find(W<Zmin);
        i=[i1;i2];
        xSPM1=xSPM;
        xSPM2=xSPM;
        xSPM.Z=W(i)';
        [i,j,k]=ind2sub(V.dim,i);
        xSPM.XYZ=[i j k]';
        xSPM.XYZmm=xSPM.M(1:3,:)*[xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
        xSPM1.Z=W(i1)';
        [i,j,k]=ind2sub(V.dim,i1);
        xSPM1.XYZ=[i j k]';
        xSPM1.XYZmm=xSPM1.M(1:3,:)*[xSPM1.XYZ; ones(1,size(xSPM1.XYZ,2))];
        xSPM2.Z=W(i2)';
        [i,j,k]=ind2sub(V.dim,i2);
        xSPM2.XYZ=[i j k]';
        xSPM2.XYZmm=xSPM2.M(1:3,:)*[xSPM2.XYZ; ones(1,size(xSPM2.XYZ,2))];
    end
end

tmp=spm_input('T1 template overlay ','+1','Yes|No');
if strcmp(tmp,'Yes')
    try
        Species=D.Atlas;
    catch
        Species=spm_input('Species ','+1','Human|Rat|Mouse');
    end
    switch Species
        case{'Human'}
            spms=fullfile(spm('dir'),'templates','T1.nii');
            spms=fullfile(spm('dir'),'canonical','single_subj_T1.nii');
            spms=fullfile(spm('dir'),'canonical','avg152T1.nii');
        case{'Rat'}
            spms=fullfile(spm_str_manip(spm('dir'),'h'),'ratlas5','template','template_T1.img');
        case{'Mouse'}
            spms=fullfile(spm_str_manip(spm('dir'),'h'),'mousatlas5','template','template_T1.img');
    end
end

%xSPM.Z=image (1*number of voxel)
%xSPM.XYZmm=location (3*number of voxel)
%M=transformation matrix
%DIM=image dimension in voxels (3*1)
spm_figure('Clear',Finter)
hReg      = spm_results_ui('SetupGUI',xSPM.M,DIM,xSPM,Finter);

spm_figure('Clear',Fgraph)
hMIPax = axes('Parent',Fgraph,'Position',[0.13 0.52 0.715 0.468],'Visible','off');
hMIPax = spm_mip_ui(abs(xSPM.Z),xSPM.XYZmm,xSPM.M,DIM,hMIPax);
% h=get(hMIPax,'Children');
% tmp=get(h(1),'Position')
% set(h(1),'Position',[tmp(1:2) 1000])
% tmp=get(h(2),'Position')
% set(h(2),'Position',[tmp(1:2) 1000])
% tmp=get(h(3),'Position')
% set(h(3),'Position',[tmp(1:2) 1000])
% tmp=get(h(4),'Position')
% set(h(4),'Position',[tmp(1:2) 1000])
% tmp=get(h(5),'Position')
% set(h(5),'Position',[tmp(1:2) 1000])
% tmp=get(h(6),'Position')
% set(h(6),'Position',[tmp(1:2) 1000])


if FlagSyn
    if ~exist('Species')
        try
            Species=D.Atlas;
        catch
            Species=spm_input('Species ','+1','Human|Rat|Mouse');
        end
    end


    %-Axis offsets for 3d MIPs:
    %=======================================================================
    %-MIP pane dimensions and Talairach origin offsets
    %-See spm_project.c for derivation
    % DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
    %-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
    % Transverse: [Po(1), Po(2)]
    % Saggital  : [Po(1), Po(3)]
    % Coronal   : [Po(4), Po(3)]
    % 4 voxel offsets in Y since using character '<' as a pointer.
    switch Species
        case{'Human'}
            DXYZ = [182 218 182];
            CXYZ = [091 127 073];
%             DXYZ = [192 218 182];
%             CXYZ = [096 132 078];
        case{'Rat','Mouse'}
            DXYZ = [189 211 175];
            CXYZ = [106 167 141];
    end
    Po(1)  =                  CXYZ(2) -2;
    Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
    Po(3)  = DXYZ(3)         -CXYZ(3) +2;
    Po(4)  = DXYZ(2)         +CXYZ(1) -2;

    axes(hMIPax)
    for i1=1:length(PairPos)
        ind=M(:,PairPos(i1));
        h=line(Po(1) + Syn.Pos(ind,2),Po(2) + Syn.Pos(ind,1));
        set(h,'Linewidth',1,'color','r');
        h=line(Po(1) + Syn.Pos(ind,2),Po(3) - Syn.Pos(ind,3));
        set(h,'Linewidth',1,'color','r');
        h=line(Po(4) + Syn.Pos(ind,1),Po(3) - Syn.Pos(ind,3));
        set(h,'Linewidth',1,'color','r');
    end
    for i1=1:length(PairNeg)
        ind=M(:,PairNeg(i1));
        h=line(Po(1) + Syn.Pos(ind,2),Po(2) + Syn.Pos(ind,1));
        set(h,'Linewidth',1);
        h=line(Po(1) + Syn.Pos(ind,2),Po(3) - Syn.Pos(ind,3));
        set(h,'Linewidth',1);
        h=line(Po(4) + Syn.Pos(ind,1),Po(3) - Syn.Pos(ind,3));
        set(h,'Linewidth',1);
    end
end

if FlagElec
    if ~exist('Species')
        try
            Species=D.Atlas;
        catch
            Species=spm_input('Species ','+1','Human|Rat|Mouse');
        end
    end
    %-Axis offsets for 3d MIPs:
    %=======================================================================
    %-MIP pane dimensions and Talairach origin offsets
    %-See spm_project.c for derivation
    % DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
    %-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
    % Transverse: [Po(1), Po(2)]
    % Saggital  : [Po(1), Po(3)]
    % Coronal   : [Po(4), Po(3)]
    % 4 voxel offsets in Y since using character '<' as a pointer.
    switch Species
%         case{'Human'}
%             DXYZ = [182 218 182];
%             CXYZ = [091 127 073];
% %             DXYZ = [192 218 182];
% %             DXYZ = [189 218 185];
% %             CXYZ = [096 128 078];
%         case{'Rat'}
%             DXYZ = [189 211 175];
%             CXYZ = [106 167 141];
        case{'Human'}
            DXYZ = [182 218 182];
            CXYZ = [091 127 073];
%             DXYZ = [192 218 182];
%             CXYZ = [096 132 078];
        case{'Rat','Mouse'}
            DXYZ = [189 211 175];
            CXYZ = [106 167 141];
    end
    Po(1)  =                  CXYZ(2) -2;
    Po(1)  =                  CXYZ(2) +3;
    Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
    Po(3)  = DXYZ(3)         -CXYZ(3) +2;
    Po(3)  = DXYZ(3)         -CXYZ(3) -3;
    Po(4)  = DXYZ(2)         +CXYZ(1) -2;

    for i0=1:size(FileElec,1)
        PosElec=load(deblank(FileElec(i0,:)));
        col=spm_input(sprintf('Color for electrode set %d ',i0),'+1','b|g|r|c|m|y|k');

        axes(hMIPax)
        for i1=1:size(PosElec,1)
            hold on
            h1=plot(Po(1) + PosElec(i1,2),Po(2) + PosElec(i1,1),'.');
            h2=plot(Po(1) + PosElec(i1,2),Po(3) - PosElec(i1,3),'.');
            h3=plot(Po(4) + PosElec(i1,1),Po(3) - PosElec(i1,3),'.');
            hold off
            colordot=col;
            if size(PosElec,1)<513
                set(h1,'MarkerSize',10,'color',colordot);
                set(h2,'MarkerSize',10,'color',colordot);
                set(h3,'MarkerSize',10,'color',colordot);
            else
                set(h1,'MarkerSize',3,'color',colordot);
                set(h2,'MarkerSize',3,'color',colordot);
                set(h3,'MarkerSize',3,'color',colordot);
            end
        end
    end
end

spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');

if exist('spms')
    if ~exist('xSPM1')
%         if FlagElec
%             xSPM.PosElec=PosElec;
%         end
        ImaGIN_Disp3D_spm_sections(xSPM,hReg,spms);
    else
%         if FlagElec
%             xSPM1.PosElec=PosElec;
%         end
        ImaGIN_Disp3D_spm_sections2(xSPM1,xSPM2,hReg,spms);
    end
else
    if ~exist('xSPM1')
%         if FlagElec
%             xSPM.PosElec=PosElec;
%         end
        ImaGIN_Disp3D_spm_sections(xSPM,hReg);
    else
%         if FlagElec
%             xSPM1.PosElec=PosElec;
%         end
        ImaGIN_Disp3D_spm_sections2(xSPM1,xSPM2,hReg);
    end
end



function ImaGIN_Disp3D_spm_sections(SPM,hReg,spms)
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,hReg)
%
% SPM  - xSPM structure containing details of excursion set
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_sections.m 112 2005-05-04 18:20:52Z john $


if nargin < 3 | isempty(spms)
	spms   = spm_select(1,'image','select image for rendering on');
end

Fgraph = spm_figure('FindWin','Graphics');
% spm_results_ui('Clear',Fgraph);
spm_orthviews('Reset');
global st
st.Space = spm_matrix([0 0 0  0 0 -pi/2])*st.Space;
spm_orthviews('Image',spms,[0.05 0.05 0.9 0.45]);
spm_orthviews MaxBB;
spm_orthviews('register',hReg);
spm_orthviews('addblobs',1,SPM.XYZ,SPM.Z,SPM.M);
spm_orthviews('Redraw');

if isfield(SPM,'PosElec')
    spm_orthviews('Reposition',SPM.PosElec(1,:)');
    colors(1)='c';
    colors(2)='m';
    for i1=1:size(SPM.PosElec,1)
        dip_h(:,i1)=add1elec(SPM.PosElec(i1,:)','.',colors(i1),st.vols{1}.ax,Fgraph,st.bb);
    end
    st.vols{1}.sdip.ax = dip_h;
end

global prevsect
prevsect = spms;





function ImaGIN_Disp3D_spm_sections2(SPM1,SPM2,hReg,spms)
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections(SPM,hReg)
%
% SPM  - xSPM structure containing details of excursion set
% hReg - handle of MIP register
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in SPM and
% VOL to create three orthogonal sections though a background image.
% Regional foci from the selected SPM are rendered on this image.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_sections.m 112 2005-05-04 18:20:52Z john $


if nargin < 4 | isempty(spms)
	spms   = spm_select(1,'image','select image for rendering on');
end

Fgraph = spm_figure('FindWin','Graphics');
% spm_results_ui('Clear',Fgraph);
spm_orthviews('Reset');
global st
st.Space = spm_matrix([0 0 0  0 0 -pi/2])*st.Space;
spm_orthviews('Image',spms,[0.05 0.05 0.9 0.45]);
spm_orthviews MaxBB;
spm_orthviews('register',hReg);
% 		c = spm_input('Colour','+1','m','Red blobs|Yellow blobs|Green blobs|Cyan blobs|Blue blobs|Magenta blobs',[1 2 3 4 5 6],1);
% 		colours = [1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1];
% 		spm_orthviews('addcolouredblobs',1,VOL.XYZ,VOL.Z,VOL.M,colours(c,:));
spm_orthviews('addcolouredblobs',1,SPM1.XYZ,SPM1.Z,SPM1.M,[1 0 0]);
spm_orthviews('Redraw');
% spm_orthviews('addcolouredblobs',1,SPM2.XYZ,SPM2.Z,SPM2.M,[1 1 0]);
spm_orthviews('addcolouredblobs',1,SPM2.XYZ,SPM2.Z,SPM2.M,[0 0 1]);
spm_orthviews('Redraw');


if isfield(SPM1,'PosElec')
    spm_orthviews('Reposition',SPM1.PosElec(1,:)');
    colors(1)='c';
    colors(2)='m';
    for i1=1:size(SPM1.PosElec,1)
        dip_h(:,i1)=add1elec(SPM1.PosElec(i1,:)','.',colors(i1),st.vols{1}.ax,Fgraph,st.bb);
    end
    st.vols{1}.sdip.ax = dip_h;
end


global prevsect
prevsect = spms;






function dh = add1elec(loc,mark,col,ax,Fig,bb)
% Plots the dipoles on the 3 views
% Then returns the handle to the plots

% loc(1,:) = loc(1,:) - bb(1,1)+1;
% loc(2,:) = loc(2,:) - bb(1,2)+1;
% loc(3,:) = loc(3,:) - bb(1,3)+1;

loc1=loc;

loc(1,:) = loc1(1,:) - bb(1,2)+1;
loc(2,:) = loc1(2,:) - bb(1,1)+1;
loc(3,:) = loc1(3,:) - bb(1,3)+1;
% +1 added to be like John's orthview code

dh = zeros(3,1);
figure(Fig)
% Transverse slice, # 1
set(Fig,'CurrentAxes',ax{1}.ax)
set(ax{1}.ax,'NextPlot','add')
% dh(1) = plot(loc(1),loc(2),[mark,col],'MarkerSize',20);
dh(1) = plot(loc(2),bb(2,2)-bb(1,2)-loc(1),[mark,col],'MarkerSize',20);
set(ax{1}.ax,'NextPlot','replace')

% Sagital slice, # 2
set(Fig,'CurrentAxes',ax{2}.ax)
set(ax{2}.ax,'NextPlot','add')
% dh(2) = plot(loc(1),loc(3),[mark,col],'MarkerSize',20);
dh(2) = plot(loc(2),loc(3),[mark,col],'MarkerSize',20);
set(ax{2}.ax,'NextPlot','replace')

% Coronal slice, # 3
set(Fig,'CurrentAxes',ax{3}.ax)
set(ax{3}.ax,'NextPlot','add')
% dh(3) = plot(bb(2,2)-bb(1,2)-loc(2),loc(3),[mark,col],'MarkerSize',20);
dh(3) = plot(loc(1),loc(3),[mark,col],'MarkerSize',20);
set(ax{3}.ax,'NextPlot','replace')

