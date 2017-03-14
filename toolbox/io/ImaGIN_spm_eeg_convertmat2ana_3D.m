function ImaGIN_spm_eeg_convertmat2ana_3D(S)

%3D for intracerebral EEG
%Olivier David

% Convert epoched EEG/ERP data from SPM- to analyze format by projecting
% onto the scalp surface
% FORMAT spm_eeg_convertmat2ana(S)
%
% S		    - optinal input struct
% (optional) fields of S:
% Fname		- matrix of EEG mat-files
% n         - size of quadratic output image (size: n x n x 1)
%_______________________________________________________________________
%
% spm_eeg_convertmat2ana converts EEG/MEG data from the SPM format to the
% scalp format. The channel data is interpolated to voxel-space using a
% spline interpolation. The electrodes' locations are specified by the
% channel template file. Each channel's data will be found in an individual
% voxel given that n is big enough. The data is written to 4-dim analyze
% images, i.e. the data of each single trial or ERP is contained in one
% image file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_convertmat2ana.m 213 2005-08-22 12:43:29Z stefan $

% [Finter, Fgraph, CmdLine] = spm('FnUIsetup', 'EEG conversion setup',0);
% 
% select matfiles to convert

try
    Fname = S.Fname;
catch
    Fname = spm_select(inf, '\.mat$', 'Select EEG mat file');
end

Nsub = size(Fname, 1);

try
    CorticalMesh = S.CorticalMesh;
catch
    str   = 'Use cortical mesh ';
    str=spm_input(str, '+1','Yes|No');
    if strcmp(str,'Yes')
        CorticalMesh = 1;
    else
        CorticalMesh = 0;
    end
end
  
if CorticalMesh
    try 
        sMRI = S.sMRI;
    catch
        sMRI = spm_select(1, 'image', 'Select normalised MRI');
    end
    try 
        SaveMNI = S.SaveMNI;
    catch
        str   = 'Save image in';
        % 	Sel   = spm_input(str, '+1', 'm', Ctype);
        % 	S.Atlas = Ctype{Sel};
        SaveMNI=spm_input(str, '+1','Native|MNI');
        if strcmp(SaveMNI,'Native')
            SaveMNI=0;
        else
            SaveMNI=1;
        end
    end
else
    try
        SizeSphere = S.SizeSphere;
    catch
        SizeSphere = spm_input('Size of local spheres [mm]', '+1', 'n', '5', 1);
    end
end

try
    SizeHorizon = S.SizeHorizon;
catch
    SizeHorizon = spm_input('Size of spatial horizon [mm]', '+1', 'n', '10', 1);
end

try
    n = S.n;
catch
    n = spm_input('Output image spatial resolution [mm]', '+1', 'n', '4', 1);
end

try
    TimeWindow = S.TimeWindow;
catch
    TimeWindow = spm_input('Time window positions [sec]', '+1', 'r');
end
if isempty(TimeWindow)
    TimeWindowFlag=1;
else
    TimeWindowFlag=0;
end

try
    TimeWindowWidth = S.TimeWindowWidth;
catch
    TimeWindowWidth = spm_input('Time window width [sec]', '+1', 'r');
end

try
    Atlas = S.Atlas;
catch
    str   = 'Select atlas';
% 	Sel   = spm_input(str, '+1', 'm', Ctype);
% 	S.Atlas = Ctype{Sel};
    Atlas=spm_input(str, '+1','Human|Rat|Mouse|PPN');
end

try
    FileOut=S.FileOut;
catch
    FileOut=[];
end


if length(n) > 1
    error('Output image spatial resolution must be scalar');
end

if length(TimeWindowWidth) > 1
    error('Time window width must be scalar');
end

try
    interpolate_bad = S.interpolate_bad;
catch
    interpolate_bad = spm_input('Interpolate bad channels or mask out?',...
        '+1', 'b', 'Interpolate|Mask out', [1,0]);
end

spm('Pointer', 'Watch'); drawnow

% Load data set into structures
clear D
for i = 1:Nsub
    D{i} = spm_eeg_load(deblank(Fname(i,:)));
end

for k = 1:Nsub
    
    if TimeWindowFlag
        TimeWindow=D{k}.time;
    end
    
    % load channel template file (contains location of channels)
%     Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
%     Ctf=load(Fchannels);
%     Ctf = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));
    Ctf=sensors(D{k},'EEG');
    
    
    if CorticalMesh

        mesh = ImaGIN_spm_eeg_inv_mesh(sMRI, 4);
        mm        = export(gifti(mesh.tess_ctx),'patch');
        mm_mni        = export(gifti(mesh.tess_mni),'patch');
        GL      = spm_mesh_smooth(mm);
        DistMean = 3.4;     %average distance between neighbours at highest resolution
        
        
        Bad=badchannels(D{k});
        Good=setdiff(setdiff(1:nchannels(D{k}),indchantype(D{k},'ECG')),Bad);
        
        Index=cell(1,nchannels(D{k})-length(indchantype(D{k},'ECG')));
        Distance=cell(1,nchannels(D{k})-length(indchantype(D{k},'ECG')));
        for i1=Good
            d=sqrt(sum((mm.vertices-ones(size(mm.vertices,1),1)*Ctf.elecpos(i1,:)).^2,2));
            Distance{i1}=min(d);
            if Distance{i1}<=SizeHorizon
                Index{i1}=find(d==min(d))';
                Distance{i1}=Distance{i1}*ones(1,length(Index{i1}));
            end
        end
        
        ok=1;
%         n=0;
        while ok
%             n=n+1
            ok=0;
            IndexConn=cell(1,length(Index));
            IndexNew=cell(1,length(Index));
            DistanceNew=cell(1,length(Index));
            %Croissance dans un volume
            for i1=Good
                for i2=1:length(Index{i1})
                    IndexConn{i1}=unique([IndexConn{i1} find(GL(Index{i1}(i2),:))]);
                end
                IndexNew{i1}=setdiff(IndexConn{i1},Index{i1});
                d=sqrt(sum((mm.vertices(IndexNew{i1},:)-ones(length(IndexNew{i1}),1)*Ctf.elecpos(i1,:)).^2,2));
                DistanceNew{i1}=d';
                DistanceNew{i1}=DistanceNew{i1}(find(d<=SizeHorizon));
                IndexNew{i1}=IndexNew{i1}(find(d<=SizeHorizon));
                if ~isempty(IndexNew{i1})
                    ok=1;
                    Index{i1}=[Index{i1} IndexNew{i1}];
                    Distance{i1}=[Distance{i1} DistanceNew{i1}];
                end
            end
        end
        Cind=Good;

    else

        [Cel, Cind, x, y, z, Index] = ImaGIN_spm_eeg_locate_channels(D{k}, n, interpolate_bad,SizeHorizon,Ctf,Atlas);
        [Cel2, Cind2, x2, y2, z2, Index2] = ImaGIN_spm_eeg_locate_channels(D{k}, n, interpolate_bad,SizeSphere,Ctf,Atlas);

    end





    if isfield(D{k},'time')
        time=D{k}.time;
    else
        time=0:1/D{k}.fsample:(D{k}.nsamples-1)/D{k}.fsample;
        time=time+D{k}.timeonset;
    end
    timewindow=TimeWindow;
    for i1=1:length(timewindow)
        [tmp,timewindow(i1)]=min(abs(time-TimeWindow(i1)));
    end
    timewindow=unique(timewindow);
%     if isfield(D{k},'time')
%     if D{k}.nsamples>1
        timewindowwidth=round(TimeWindowWidth*D{k}.fsample/2);
%     else
%         timewindowwidth=0;
%     end
    
    switch Atlas
        case{'Human'}
            bb = [[-78 -112 -50];[78 76 85]];
            tmp=spm('Defaults','EEG');
            bb=tmp.normalise.write.bb;
            V = fullfile(spm('dir'), 'toolbox', 'OldNorm', 'T1.nii');
            V=spm_vol(V);
        case{'PPN'}
            bb = [[-8 -5 -20];[8 6 2]];     %Brainstem full
            V = '/Users/odavid/Documents/Data/Goetz/IRM/MRI_PPN_Small2.img';
            V=spm_vol(V);
        case{'Rat'}
            bb = [[-80 -156 -120];[80 60 10]];
            V = fullfile(spm('dir'),'atlas8','rat','template','template_T1.img');
            V=spm_vol(V);
        case{'Mouse'}
            bb = [[-48 -94 -70];[48 72 0]];
            V = fullfile(spm('dir'),'atlas8','mouse','template','template_T1.img');
            V=spm_vol(V);
    end
    n1=length(bb(1,1):n:bb(2,1));
    n2=length(bb(1,2):n:bb(2,2));
    n3=length(bb(1,3):n:bb(2,3));

    % generate data directory into which converted data goes
    [P, F] = fileparts(spm_str_manip(Fname(k, :), 'r'));
    if ~isempty(P)
        [m, sta] = mkdir(P, spm_str_manip(Fname(k, :), 'tr'));
    else
        mkdir(spm_str_manip(Fname(k, :), 'tr'));
    end
    cd(fullfile(P, F));

    %Check if it is synchrony
    if isempty(strfind(D{k}.fname,'2int'))
        FlagSyn=0;
        d = (D{k}(Cind, :,:));
    else
        FlagSyn=1;
%         Nchannels=length(D{k}.channels.eeg)-length(D{k}.channels.Bad);
        Nchannels=(1+sqrt(1+8*D{k}.nchannels))/2;
        M=ImaGIN_ConnectivityMatrix(Nchannels);
        tmpd=(D{k}(:, :,:));
        d=zeros(Nchannels,D{k}.nsamples);
        for i1=1:Nchannels
            [tmp1,tmp2]=find(M==i1);
            d(i1,:)=mean(tmpd(tmp2,:));
        end
        if sign(min(tmpd(:)))==sign(max(tmpd(:)))
            d1=zeros(Nchannels,D{k}.nsamples);
            d2=zeros(Nchannels,D{k}.nsamples);
            for i1=1:Nchannels
                [tmp1,tmp2]=find(M==i1);
                for i2=1:size(tmpd,2)
                    tmp1=find(tmpd(tmp2,i2)>=0);
                    if ~isempty(tmp1)
                        d1(i1,:)=mean(tmpd(tmp2(tmp1),i2),1);
                    end
                    tmp1=find(tmpd(tmp2,i2)<0);
                    if ~isempty(tmp1)
                        d2(i1,:)=mean(tmpd(tmp2(tmp1),i2),1);
                    end
                end
            end
        end
    end

%     D{k}.Nsamples=2;
    
%     fname = [F '.img'];
%     dat = file_array(fname,[n1 n2 n3 D{k}.Nsamples],'FLOAT32');
%     N = nifti;
%     N.dat = dat;
%     N.mat = eye(4);
%     N.mat_intent = 'Aligned';
%     create(N);
%     for j = 1 : D{k}.Nsamples % time bins
%         di = zeros(n1, n2, n3);
%         di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), d(:, j),x,y,z, 'linear');
%         N.dat(:,:,:,j) = di;
%         disp(sprintf('Subject %d, sample %d', k, j))
%     end

    tmp=round(1000*max(abs(time(timewindow))));
    for j = timewindow % time bins
        J=round(1000*time(j));
        if tmp<1e1
            V.fname = sprintf('sample_%d.img',J);
        elseif tmp<1e2
            V.fname = sprintf('sample_%0.2d.img',J);
        elseif tmp<1e3
            V.fname = sprintf('sample_%0.3d.img',J);
        elseif tmp<1e4
            V.fname = sprintf('sample_%0.4d.img',J);
        elseif tmp<1e5
            V.fname = sprintf('sample_%0.5d.img',J);
        elseif tmp<1e6
            V.fname = sprintf('sample_%0.6d.img',J);
        else
            V.fname = sprintf('sample_%d.img',J);            
        end            
%         if D{k}.Nsamples<10
%             V.fname = sprintf('sample_%d.img',J);
%         elseif D{k}.Nsamples<100
%             V.fname = sprintf('sample_%0.2d.img',J);
%         elseif D{k}.Nsamples<1000
%             V.fname = sprintf('sample_%0.3d.img',J);
%         elseif D{k}.Nsamples<10000
%             V.fname = sprintf('sample_%d.img',J);
%         end            
        
        win=j+[-timewindowwidth:timewindowwidth];
        win=win(find(win>=1&win<=D{k}.nsamples));
        tmpd=mean(d(:,win),2);
        
        P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
        V.mat=spm_matrix(P);
        V.dim=[n1 n2 n3];
        if isfield(S,'dt')
            V.dt=S.dt;
        else
            V.dt=[64 0];    %float 64
            V.dt=[16 0];    %float 32
            V.dt=[16 0];    %float 32
        end
        
        if CorticalMesh
            
            EMap=zeros(length(GL),1);
            EMapDist=zeros(length(GL),1);
            for i1=1:length(Cind)
                if isnan(tmpd(i1))
                    map=EMapDist(Index{Cind(i1)});
                    mapZero=find(map==0);
                    EMap(Index{Cind(i1)}(mapZero))=NaN;
                else
                    map=EMap(Index{Cind(i1)});
                    mapNoNaN=find(~isnan(map));
                    mapNaN=find(isnan(map));
                    EMap(Index{Cind(i1)}(mapNoNaN))=EMap(Index{Cind(i1)}(mapNoNaN))+tmpd(i1)*(SizeHorizon-Distance{Cind(i1)}(mapNoNaN))';
                    EMapDist(Index{Cind(i1)}(mapNoNaN))=EMapDist(Index{Cind(i1)}(mapNoNaN))+SizeHorizon-Distance{Cind(i1)}(mapNoNaN)';
                    EMap(Index{Cind(i1)}(mapNaN))=tmpd(i1)*(SizeHorizon-Distance{Cind(i1)}(mapNaN))';
                    EMapDist(Index{Cind(i1)}(mapNaN))=SizeHorizon-Distance{Cind(i1)}(mapNaN)';
                end
            end
            EMap=EMap./EMapDist;
            
            if SaveMNI
                di = ImaGIN_spm_mesh_to_grid(mm_mni, V, EMap);
            else
                di = ImaGIN_spm_mesh_to_grid(mm, V, EMap);
            end
%             Maxdi=max(di(:));
%             if abs(Maxdi-1)<1e-6
%                 V.dt=[4 0];
%             end
            Maskout=find(isnan(di));
            di(Maskout)=-Inf;
            di=spm_dilate(di);
%             di2=spm_dilate(abs(di));
%             di3=sign(di);di3(di3==0)=1;
%             di3=spm_dilate(di3);
%             di=di2.*di3;
%             di=spm_dilate(di);
            di(di==-Inf)=NaN;
% %             di=(Maxdi/max(di(:)))*di;
% %             
% %         spm_write_vol(V,di);
% % %         V=spm_vol(V.fname);
% % %         spm_smooth(V,V.fname,2);
% %             %smooth image of latency
% %             clear matlabbatch
% %             matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
% %             matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
% %             matlabbatch{1}.spm.spatial.smooth.dtype = 0;
% %             matlabbatch{1}.spm.spatial.smooth.im = 1;
% %             matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% %             spm('defaults', 'EEG');
% %             spm_jobman('run', matlabbatch);
% %             movefile(['s' V.fname],V.fname)
% %             movefile(['s' V.fname(1:end-3) 'hdr'],[V.fname(1:end-3) 'hdr'])
% %             
% %             
% %             di=spm_vol(V.fname);
% %             di=spm_read_vols(di);
% %             di=(Maxdi/max(di(:)))*di;
%             %erosion
%             Mask=zeros(size(di));
%             Mask(~isnan(di))=1;
%             Mask=spm_erode(Mask);
%             Mask(Mask==0)=NaN;
%             di=di.*Mask;
%             di=(Maxdi/max(di(:)))*di;

            V=spm_write_vol(V,di);

%             clear matlabbatch
%             matlabbatch{1}.spm.spatial.smooth.data = {[V.fname ',1']};
%             matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
%             matlabbatch{1}.spm.spatial.smooth.dtype = 0;
%             matlabbatch{1}.spm.spatial.smooth.im = 1;
%             matlabbatch{1}.spm.spatial.smooth.prefix = 's';
%             spm('defaults', 'EEG');
%             spm_jobman('run', matlabbatch);
%             movefile(['s' V.fname],V.fname)
%             movefile(['s' V.fname(1:end-3) 'hdr'],[V.fname(1:end-3) 'hdr'])
% 
%             di=spm_vol(V.fname);
%             di=spm_read_vols(di);
%             di=(Maxdi/max(di(:)))*di;
%             V=spm_write_vol(V,di);
% 
%             %             spm_smooth(di, di, 0.5);
% %             di = di.*(abs(di) > Maxdi*exp(-8));            
% %             di(abs(di) <= Maxdi*exp(-8))=NaN;
% %             di=(Maxdi/max(di(:)))*di;
% %             %erosion
% %             Mask=zeros(size(di));
% %             Mask(abs(di)>Maxdi*exp(-8))=1;
% %             Mask=spm_erode(Mask);
% %             Mask=spm_erode(Mask);
% %             di=di.*Mask;
%             
%             
% %             di=spm_erode(di);
% %             di=spm_erode(di);
            

        else
            
            di = NaN*zeros(n2,n1,n3);
            %         di = zeros(n2,n1,n3);
            di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), tmpd,x,y,z, 'linear');
            %         di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), tmpd,x,y,z, 'linear',{'Qt','Qbb','Qc','Qz'});
            %         di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), tmpd,x,y,z, 'nearest');
            di2 = zeros(n2,n1,n3);
            di2 = NaN*zeros(n2,n1,n3);
            di2(Index2) = griddata3(Cel2(:,1), Cel2(:,2), Cel2(:,3), tmpd,x2,y2,z2, 'nearest');
            Index3=intersect(find(isnan(di)),Index2);
            di(Index3)=di2(Index3);
            %         di(Index)=1;
            di=permute(di,[2 1 3]);
            
        V=spm_write_vol(V,di);
%         V=spm_vol(V.fname);
%         spm_smooth(V,V.fname,2);
        end
        
        
        
        
        
        if exist('d1')
            tmpd=mean(d1(:,win),2);
            di = NaN*zeros(n2,n1,n3);
            di = zeros(n2,n1,n3);
            di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), tmpd,x,y,z, 'nearest');
            %         di(Index)=1;
            di=permute(di,[2 1 3]);
            P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
            Vtmp=V;
            Vtmp.mat=spm_matrix(P);
            Vtmp.dim=[n1 n2 n3];
            Vtmp.dt=[64 0];
            Vtmp.dt=[16 0];
            Vtmp.fname=['pos_' V.fname];
            Vtmp=spm_write_vol(Vtmp,di);
%             Vtmp=spm_vol(Vtmp.fname);
%             spm_smooth(Vtmp,Vtmp.fname,2);
            tmpd=mean(d2(:,win),2);
            di = NaN*zeros(n2,n1,n3);
            di = zeros(n2,n1,n3);
            di(Index) = griddata3(Cel(:,1), Cel(:,2), Cel(:,3), tmpd,x,y,z, 'nearest');
            %         di(Index)=1;
            di=permute(di,[2 1 3]);
            P=[bb(1,1),bb(1,2),bb(1,3),0,0,0,n,n,n];
            Vtmp=V;
            Vtmp.mat=spm_matrix(P);
            Vtmp.dim=[n1 n2 n3];
            Vtmp.dt=[64 0];
            Vtmp.dt=[16 0];
            Vtmp.fname=['neg_' V.fname];
            Vtmp=spm_write_vol(Vtmp,di);
%             Vtmp=spm_vol(Vtmp.fname);
%             spm_smooth(Vtmp,Vtmp.fname,2);
        end
        if FlagSyn
            dsyn=(D{k}(:, :,:));
            S=mean(dsyn(:,win),2);
            Pos=Cel;
            file=fullfile(FileOut,[spm_str_manip(V.fname,'r') '.syn']);
            save(file,'S','Pos')
            Vtmp=V;
            Vtmp.dt=[64 0];
            Vtmp.dt=[16 0];
            Vtmp.mat=spm_matrix(P);
            Vtmp.dim=[size(S,1) 1 1];
            Vtmp.fname=['syn_' V.fname];
            Vtmp=spm_write_vol(Vtmp,S);
        end
        disp(sprintf('Subject %d, time %d', k, J))
%         disp(sprintf('Subject %d, sample %d', k, j))
    end

    cd ..
end

spm('Pointer', 'Arrow');
