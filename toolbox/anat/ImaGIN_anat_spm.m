function ImaGIN_anat_spm(Patient)


for I = 1:length(Patient)
    
    matlabbatch = {};
    
    %% Coregister
    switch Patient{I}.MRI.ref
        case 'pre'
            
            %Initial translation according to centroids
            Vref=spm_vol(Patient{I}.MRI.pre);
            [Iref,XYZref]=spm_read_vols(Vref);
            Iindex=find(Iref>max(Iref(:))/6);
            Zindex=find(max(XYZref(3,:))-XYZref(3,:)<200);
            index=intersect(Iindex,Zindex);
            CentroidRef=mean(XYZref(:,index),2);
            if isfield(Patient{I}.MRI,'per')
                V2=spm_vol(Patient{I}.MRI.per);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    V2=spm_vol(Patient{I}.CT.per);
                    [I2,XYZ2]=spm_read_vols(V2);
                    Iindex=find(I2>max(I2(:))/6);
                    Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                    index=intersect(Iindex,Zindex);
                    Centroid2=mean(XYZ2(:,index),2);
                    %apply translation
                    B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                    M = spm_matrix(B);
                    Mat = spm_get_space(V2.fname);
                    spm_get_space(V2.fname,M*Mat);
                end
            end
            if isfield(Patient{I}.MRI,'post')
                V2=spm_vol(Patient{I}.MRI.post);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            
            if isfield(Patient{I}.MRI,'per')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.per};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    iReg = length(matlabbatch) + 1;
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.CT.per};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                end
            end
            if isfield(Patient{I}.MRI,'post')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            
        case 'per'
            
            %Initial translation according to centroids
            Vref=spm_vol(Patient{I}.MRI.per);
            [Iref,XYZref]=spm_read_vols(Vref);
            Iindex=find(Iref>max(Iref(:))/6);
            Zindex=find(max(XYZref(3,:))-XYZref(3,:)<200);
            index=intersect(Iindex,Zindex);
            CentroidRef=mean(XYZref(:,index),2);
            if isfield(Patient{I}.MRI,'pre')
                V2=spm_vol(Patient{I}.MRI.pre);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            if isfield(Patient{I}.MRI,'post')
                V2=spm_vol(Patient{I}.MRI.post);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            
            if isfield(Patient{I}.MRI,'pre')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.per};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I}.MRI,'post')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.per};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            
        case 'post'

            %Initial translation according to centroids
            Vref=spm_vol(Patient{I}.MRI.post);
            [Iref,XYZref]=spm_read_vols(Vref);
            Iindex=find(Iref>max(Iref(:))/6);
            Zindex=find(max(XYZref(3,:))-XYZref(3,:)<200);
            index=intersect(Iindex,Zindex);
            CentroidRef=mean(XYZref(:,index),2);
            if isfield(Patient{I}.MRI,'pre')
                V2=spm_vol(Patient{I}.MRI.pre);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            if isfield(Patient{I}.MRI,'per')
                V2=spm_vol(Patient{I}.MRI.per);
                [I2,XYZ2]=spm_read_vols(V2);
                Iindex=find(I2>max(I2(:))/6);
                Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                index=intersect(Iindex,Zindex);
                Centroid2=mean(XYZ2(:,index),2);
                %apply translation
                B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                M = spm_matrix(B);
                Mat = spm_get_space(V2.fname);
                spm_get_space(V2.fname,M*Mat);
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    V2=spm_vol(Patient{I}.CT.per);
                    [I2,XYZ2]=spm_read_vols(V2);
                    Iindex=find(I2>max(I2(:))/6);
                    Zindex=find(max(XYZ2(3,:))-XYZ2(3,:)<200);
                    index=intersect(Iindex,Zindex);
                    Centroid2=mean(XYZ2(:,index),2);
                    %apply translation
                    B=[CentroidRef'-Centroid2' 0 0 0 1 1 1 0 0 0];
                    M = spm_matrix(B);
                    Mat = spm_get_space(V2.fname);
                    spm_get_space(V2.fname,M*Mat);
                end
            end

            if isfield(Patient{I}.MRI,'pre')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I}.MRI,'per')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.per};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    iReg = length(matlabbatch) + 1;
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.post};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.CT.per};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                end
            end
    end

    
    %% Segmentation
    if isfield(Patient{I}.MRI,'pre')
        iSeg = length(matlabbatch) + 1;
        matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{iSeg}.spm.spatial.preproc.channel.vols = {Patient{I}.MRI.pre};
        ngaus  = [1 1 2 3 4 2];
        native = [1 1 1 0 0 0];
        for c = 1:6 % tissue class c
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
                fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
        end
        matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];
    end
    if isfield(Patient{I}.MRI,'per')
        iSeg = length(matlabbatch) + 1;
        matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{iSeg}.spm.spatial.preproc.channel.vols = {Patient{I}.MRI.per};
        ngaus  = [1 1 2 3 4 2];
        native = [1 1 1 0 0 0];
        for c = 1:6 % tissue class c
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
                fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
        end
        matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];
    end
    if isfield(Patient{I}.MRI,'post')
        iSeg = length(matlabbatch) + 1;
        matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{iSeg}.spm.spatial.preproc.channel.vols = {Patient{I}.MRI.post};
        ngaus  = [1 1 2 3 4 2];
        native = [1 1 1 0 0 0];
        for c = 1:6 % tissue class c
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).tpm = {
                fullfile(spm('dir'), 'tpm', sprintf('TPM.nii,%d', c))};
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
            matlabbatch{iSeg}.spm.spatial.preproc.tissue(c).warped = [0 0];
        end
        matlabbatch{iSeg}.spm.spatial.preproc.warp.write = [0 1];
    end
    
    %% Image Calculator - Create brain image (skull-stripped bias corrected)
    if isfield(Patient{I}.MRI,'pre')
        iIC = length(matlabbatch) + 1;
        tmp1=spm_str_manip(Patient{I}.MRI.pre,'h');
        tmp2=spm_str_manip(Patient{I}.MRI.pre,'rt');
        matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPre';
        matlabbatch{iIC}.spm.util.imcalc.input = {
            fullfile(tmp1,['c1' tmp2 '.nii'])
            fullfile(tmp1,['c2' tmp2 '.nii'])
            fullfile(tmp1,['c3' tmp2 '.nii'])
            fullfile(tmp1,['m' tmp2 '.nii'])
            };
        matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
        matlabbatch{iIC}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
        if strcmp(Patient{I}.MRI.ref,'pre')
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    iIC = length(matlabbatch) + 1;
                    tmp5=spm_str_manip(Patient{I}.CT.per,'h');
                    tmp6=spm_str_manip(Patient{I}.CT.per,'t');
                    matlabbatch{iIC}.spm.util.imcalc.input = {
                        fullfile(tmp5,tmp6)
                        fullfile(tmp1,['c1' tmp2 '.nii'])
                        fullfile(tmp1,['c2' tmp2 '.nii'])
                        };
                    matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPerCT';
                    matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
                    matlabbatch{iIC}.spm.util.imcalc.expression = '(i2 + i3) .* i1';
                end
            end
        end
    end
    if isfield(Patient{I}.MRI,'per')
        iIC = length(matlabbatch) + 1;
        tmp1=spm_str_manip(Patient{I}.MRI.per,'h');
        tmp2=spm_str_manip(Patient{I}.MRI.per,'rt');
        matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPer';
        matlabbatch{iIC}.spm.util.imcalc.input = {
            fullfile(tmp1,['c1' tmp2 '.nii'])
            fullfile(tmp1,['c2' tmp2 '.nii'])
            fullfile(tmp1,['c3' tmp2 '.nii'])
            fullfile(tmp1,['m' tmp2 '.nii'])
            };
        matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
        matlabbatch{iIC}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
    end
    if isfield(Patient{I}.MRI,'post')
        iIC = length(matlabbatch) + 1;
        tmp1=spm_str_manip(Patient{I}.MRI.post,'h');
        tmp2=spm_str_manip(Patient{I}.MRI.post,'rt');
        matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPost';
        matlabbatch{iIC}.spm.util.imcalc.input = {
            fullfile(tmp1,['c1' tmp2 '.nii'])
            fullfile(tmp1,['c2' tmp2 '.nii'])
            fullfile(tmp1,['c3' tmp2 '.nii'])
            fullfile(tmp1,['m' tmp2 '.nii'])
            };
        matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
        matlabbatch{iIC}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
        if strcmp(Patient{I}.MRI.ref,'post')
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'per')
                    iIC = length(matlabbatch) + 1;
                    tmp5=spm_str_manip(Patient{I}.CT.per,'h');
                    tmp6=spm_str_manip(Patient{I}.CT.per,'t');
                    matlabbatch{iIC}.spm.util.imcalc.input = {
                        fullfile(tmp5,tmp6)
                        fullfile(tmp1,['c1' tmp2 '.nii'])
                        fullfile(tmp1,['c2' tmp2 '.nii'])
                        };
                    matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPerCT';
                    matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
                    matlabbatch{iIC}.spm.util.imcalc.expression = '(i2 + i3) .* i1';
                end
            end
        end
    end
    
    %% Normalise Brain image
    iNB = length(matlabbatch) + 1;
    switch Patient{I}.MRI.ref
        case 'pre'
            tmp1=spm_str_manip(Patient{I}.MRI.pre,'h');
            tmp2=spm_str_manip(Patient{I}.MRI.pre,'rt');
        case 'per'
            tmp1=spm_str_manip(Patient{I}.MRI.per,'h');
            tmp2=spm_str_manip(Patient{I}.MRI.per,'rt');
        case 'post'
            tmp1=spm_str_manip(Patient{I}.MRI.post,'h');
            tmp2=spm_str_manip(Patient{I}.MRI.post,'rt');
    end
    matlabbatch{iNB}.spm.spatial.normalise.write.subj.def = {fullfile(tmp1,['y_' tmp2 '.nii'])};
    matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample = {};
    if isfield(Patient{I}.MRI,'pre')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPre.nii');
    end
    if isfield(Patient{I}.MRI,'per')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPer.nii');
    end
    if isfield(Patient{I},'CT')
        if isfield(Patient{I}.CT,'per')
            matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPerCT.nii');
        end
    end
    if isfield(Patient{I}.MRI,'post')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPost.nii');
    end
    matlabbatch{iNB}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    
    
    spm_jobman('initcfg');
    % spm_jobman('interactive', matlabbatch)
    spm_jobman('run',matlabbatch)
    
    
    if strcmp(Patient{I}.MRI.ref,'pre') && isfield(Patient{I},'CT')
        V1=spm_vol(fullfile(Patient{I}.MRI.out,'wBrainPre.nii'));
        I1=spm_read_vols(V1);
        V2=spm_vol(fullfile(Patient{I}.MRI.out,'wBrainPerCT.nii'));
        I2=spm_read_vols(V2);
        I2=I2./max(I2(:));
        Mask=find(I2>0.3);
        V3=V1;
        V3.fname=fullfile(Patient{I}.MRI.out,'wBrainPreCT.nii');
        V3.dt=[8 0];
        I3=I1;
        I3(Mask)=2*max(I1(:));%*I2(Mask);
        spm_write_vol(V3,I3);
    end
    
    
    
end
