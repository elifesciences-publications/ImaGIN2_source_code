function ImaGIN_anat_spm(Patient)
% IMAGIN_ANAT_SPM Registration, segmentation and normalization of MRI and CT scans for SEEG/ECOG implantations.
%
% USAGE: ImaGIN_anat_spm(Patient)
%
% INPUT: 
%    - Patient{}:  Cell-array of strctures, each one representing a patient, with the following optional fields
%      |- MRI.pre    : Path to the pre-implantation MRI scan  (before any surgery)
%      |- MRI.post   : Path to the post-implantation MRI scan (where the SEEG/ECOG contacts are visible)
%      |- CT.post    : Path to the post-implantation CT scan  (where the SEEG/ECOG contacts are visible)
%      |- MRI.postop : Path to the post-surgery MRI scan      (typically after a tissue resection)
%      |- MRI.ref    : String {'pre','post','postop'}, type of the image to use as the reference for the coordinates of all the images
%      |- MRI.out    : Path to the output folder, where the normalized volumes are saved by SPM
%
% OUTPUT:  Files saved in the output folder Patient{i}.MRI.out
%    - BrainPre.nii     : Registered pre-implantation MRI
%    - BrainPost.nii    : Registered post-implantation MRI
%    - BrainPostCT.nii  : Registered post-implantation CT
%    - BrainPostOp.nii  : Registered post-surgery MRI
%    - w*.nii           : Same as above, but normalized in MNI space
%    - y_*.nii          : Deformation fields from subject space to MNI space
%    - *.surf.gii       : Surfaces reconstructed by SPM in MNI space (cortex surface, inner skull, outer skull, scalp)

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
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'post')
                    V2=spm_vol(Patient{I}.CT.post);
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
            if isfield(Patient{I}.MRI,'postop')
                V2=spm_vol(Patient{I}.MRI.postop);
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
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'post')
                    iReg = length(matlabbatch) + 1;
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.CT.post};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                end
            end
            if isfield(Patient{I}.MRI,'postop')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.postop};
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
            if isfield(Patient{I}.MRI,'postop')
                V2=spm_vol(Patient{I}.MRI.postop);
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
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I}.MRI,'postop')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.postop};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            
        case 'postop'

            %Initial translation according to centroids
            Vref=spm_vol(Patient{I}.MRI.postop);
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
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'post')
                    V2=spm_vol(Patient{I}.CT.post);
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
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.postop};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.pre};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I}.MRI,'post')
                iReg = length(matlabbatch) + 1;
                matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.postop};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.MRI.post};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.other = {''};
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            end
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'post')
                    iReg = length(matlabbatch) + 1;
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.ref = {Patient{I}.MRI.postop};
                    matlabbatch{iReg}.spm.spatial.coreg.estimate.source = {Patient{I}.CT.post};
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
    if isfield(Patient{I}.MRI,'postop')
        iSeg = length(matlabbatch) + 1;
        matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
        matlabbatch{iSeg}.spm.spatial.preproc.channel.vols = {Patient{I}.MRI.postop};
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
                if isfield(Patient{I}.CT,'post')
                    iIC = length(matlabbatch) + 1;
                    tmp5=spm_str_manip(Patient{I}.CT.post,'h');
                    tmp6=spm_str_manip(Patient{I}.CT.post,'t');
                    matlabbatch{iIC}.spm.util.imcalc.input = {
                        fullfile(tmp5,tmp6)
                        fullfile(tmp1,['c1' tmp2 '.nii'])
                        fullfile(tmp1,['c2' tmp2 '.nii'])
                        };
                    matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPostCT';
                    matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
                    matlabbatch{iIC}.spm.util.imcalc.expression = '(i2 + i3) .* i1';
                end
            end
        end
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
    end
    if isfield(Patient{I}.MRI,'postop')
        iIC = length(matlabbatch) + 1;
        tmp1=spm_str_manip(Patient{I}.MRI.postop,'h');
        tmp2=spm_str_manip(Patient{I}.MRI.postop,'rt');
        matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPostOp';
        matlabbatch{iIC}.spm.util.imcalc.input = {
            fullfile(tmp1,['c1' tmp2 '.nii'])
            fullfile(tmp1,['c2' tmp2 '.nii'])
            fullfile(tmp1,['c3' tmp2 '.nii'])
            fullfile(tmp1,['m' tmp2 '.nii'])
            };
        matlabbatch{iIC}.spm.util.imcalc.outdir = {Patient{I}.MRI.out};
        matlabbatch{iIC}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
        if strcmp(Patient{I}.MRI.ref,'postop')
            if isfield(Patient{I},'CT')
                if isfield(Patient{I}.CT,'post')
                    iIC = length(matlabbatch) + 1;
                    tmp5=spm_str_manip(Patient{I}.CT.post,'h');
                    tmp6=spm_str_manip(Patient{I}.CT.post,'t');
                    matlabbatch{iIC}.spm.util.imcalc.input = {
                        fullfile(tmp5,tmp6)
                        fullfile(tmp1,['c1' tmp2 '.nii'])
                        fullfile(tmp1,['c2' tmp2 '.nii'])
                        };
                    matlabbatch{iIC}.spm.util.imcalc.output = 'BrainPostCT';
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
        case 'post'
            tmp1=spm_str_manip(Patient{I}.MRI.post,'h');
            tmp2=spm_str_manip(Patient{I}.MRI.post,'rt');
        case 'postop'
            tmp1=spm_str_manip(Patient{I}.MRI.postop,'h');
            tmp2=spm_str_manip(Patient{I}.MRI.postop,'rt');
    end
    matlabbatch{iNB}.spm.spatial.normalise.write.subj.def = {fullfile(tmp1,['y_' tmp2 '.nii'])};
    matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample = {};
    if isfield(Patient{I}.MRI,'pre')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPre.nii');
    end
    if isfield(Patient{I}.MRI,'post')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPost.nii');
    end
    if isfield(Patient{I},'CT')
        if isfield(Patient{I}.CT,'post')
            matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPostCT.nii');
        end
    end
    if isfield(Patient{I}.MRI,'postop')
        matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPostOp.nii');
    end
    matlabbatch{iNB}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    
    % Save SPM batch
    save(fullfile(Patient{I}.MRI.out, 'ImaGIN_spm_batch.mat'), 'matlabbatch');
    % Run SPM batch
    spm_jobman('initcfg');
    % spm_jobman('interactive', matlabbatch)
    spm_jobman('run',matlabbatch)
    
    
    if strcmp(Patient{I}.MRI.ref,'pre') && isfield(Patient{I},'CT')
        V1=spm_vol(fullfile(Patient{I}.MRI.out,'wBrainPre.nii'));
        I1=spm_read_vols(V1);
        V2=spm_vol(fullfile(Patient{I}.MRI.out,'wBrainPostCT.nii'));
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


