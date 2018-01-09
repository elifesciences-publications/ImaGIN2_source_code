function Patient = ImaGIN_AnatSPM(Patient, isNormalize)
% IMAGIN_ANATSPM Registration, segmentation and normalization of MRI and CT scans for SEEG/ECOG implantations.
%
% USAGE: ImaGIN_AnatSPM(Patient, isNormalize)
%
% INPUT: 
%    - Patient{}:  Cell-array of strctures, each one representing a patient, with the following optional fields
%      |- MRI.pre        : Path to the pre-implantation MRI scan  (before any surgery)
%      |- MRI.pre_cortex : Path to the cortex surface corresponding to the "pre" MRI (in .gii format only)
%      |                   Set it to 'canonical' to compute the SPM canonical mesh
%      |- MRI.post       : Path to the post-implantation MRI scan (where the SEEG/ECOG contacts are visible)
%      |- CT.post        : Path to the post-implantation CT scan  (where the SEEG/ECOG contacts are visible)
%      |- MRI.postop     : Path to the post-surgery MRI scan      (typically after a tissue resection)
%      |- MRI.ref        : String {'pre','post','postop'}, type of the image to use as the reference for the coordinates of all the images
%      |- MRI.out        : Path to the output folder, where the normalized volumes are saved by SPM
%    - isNormalize  : If 1, coregistration pre/post/postop + skull stripping + MNI normalization
%                     If 0, coregistration pre/post/postop only
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
%          Francois Tadel, 2017

% If there is only one patient in input
if isstruct(Patient)
    Patient = {Patient};
end


for I = 1:length(Patient)
    
    matlabbatch = {};
    
    
    %% ===== COREGISTRATION =====
    % Reference for the registration
    switch Patient{I}.MRI.ref
        case 'pre'
            % Initial translation according to centroids
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
            % Register all volumes
            if isfield(Patient{I}.MRI,'post')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.pre, Patient{I}.MRI.post, isNormalize, Patient{I}.MRI.out);
            end
            if isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.pre, Patient{I}.CT.post, isNormalize, Patient{I}.MRI.out);
            end
            if isfield(Patient{I}.MRI,'postop')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.pre, Patient{I}.MRI.postop, isNormalize, Patient{I}.MRI.out);
            end
            
        case 'post'
            % Initial translation according to centroids
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
            % Register all volumes
            if isfield(Patient{I}.MRI,'pre')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.post, Patient{I}.MRI.pre, isNormalize, Patient{I}.MRI.out);
            end
            if isfield(Patient{I}.MRI,'postop')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.post, Patient{I}.MRI.postop, isNormalize, Patient{I}.MRI.out);
            end
            
        case 'postop'
            % Initial translation according to centroids
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
            % Register all volumes
            if isfield(Patient{I}.MRI,'pre')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.postop, Patient{I}.MRI.pre, isNormalize, Patient{I}.MRI.out);
            end
            if isfield(Patient{I}.MRI,'post')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.postop, Patient{I}.MRI.post, isNormalize, Patient{I}.MRI.out);
            end
            if isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
                matlabbatch = AddCoreg(matlabbatch, Patient{I}.MRI.postop, Patient{I}.CT.post, isNormalize, Patient{I}.MRI.out);
            end
    end

    % Skull stripping + MNI normalization
    if isNormalize
       
        %% ===== SEGMENTATION =====
        if isfield(Patient{I}.MRI,'pre')
            matlabbatch = AddSegment(matlabbatch, Patient{I}.MRI.pre);
        end
        if isfield(Patient{I}.MRI,'post')
            matlabbatch = AddSegment(matlabbatch, Patient{I}.MRI.post);
        end
        if isfield(Patient{I}.MRI,'postop')
            matlabbatch = AddSegment(matlabbatch, Patient{I}.MRI.postop);
        end
    
        %% ===== SKULL STRIPPING ===== 
        % Create brain image (skull-stripped bias corrected)
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
            if strcmp(Patient{I}.MRI.ref,'pre') && isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
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
            if strcmp(Patient{I}.MRI.ref,'postop') && isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
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
    
        %% ===== MNI NORMALIZATION =====
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
            Patient{I}.MRI.reg_pre = fullfile(Patient{I}.MRI.out,'wBrainPre.nii');
        end
        if isfield(Patient{I}.MRI,'post')
            matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPost.nii');
            Patient{I}.MRI.reg_post = fullfile(Patient{I}.MRI.out,'wBrainPost.nii');
        end
        if isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
            matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPostCT.nii');
            Patient{I}.MRI.reg_postct = fullfile(Patient{I}.MRI.out,'wBrainPostCT.nii');
        end
        if isfield(Patient{I}.MRI,'postop')
            matlabbatch{iNB}.spm.spatial.normalise.write.subj.resample{end+1,1} = fullfile(Patient{I}.MRI.out,'BrainPostOp.nii');
            Patient{I}.MRI.reg_postop = fullfile(Patient{I}.MRI.out,'wBrainPostOp.nii');
        end
        matlabbatch{iNB}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    else
        if isfield(Patient{I}.MRI,'pre')
            Patient{I}.MRI.reg_pre = fullfile(Patient{I}.MRI.out,'BrainPre.nii');
        end
        if isfield(Patient{I}.MRI,'post')
            Patient{I}.MRI.reg_post = fullfile(Patient{I}.MRI.out,'BrainPost.nii');
        end
        if isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
            Patient{I}.MRI.reg_postct = fullfile(Patient{I}.MRI.out,'BrainPostCT.nii');
        end
        if isfield(Patient{I}.MRI,'postop')
            Patient{I}.MRI.reg_postop = fullfile(Patient{I}.MRI.out,'BrainPostOp.nii');
        end
    end

    
    %% ===== RUN BATCH =====
    % Save SPM batch
    save(fullfile(Patient{I}.MRI.out, 'ImaGIN_spm_batch.mat'), 'matlabbatch');
    % Run SPM batch
    spm_jobman('initcfg');
    % spm_jobman('interactive', matlabbatch)
    spm_jobman('run',matlabbatch)
    
    
    %% ===== NO NORMALIZATION =====
    % When not normalizing in MNI space: Some files must be copied/moved
    if ~isNormalize
        switch Patient{I}.MRI.ref
            case 'pre'
                % Copy original PRE => BrainPre.nii
                Patient{I}.MRI.reg_pre = fullfile(Patient{I}.MRI.out, 'BrainPre.nii');
                copyfile(Patient{I}.MRI.pre, Patient{I}.MRI.reg_pre);
                % Move registered POST => BrainPost.nii
                if isfield(Patient{I}.MRI,'post')
                    [fPath, fBase, fExt] = fileparts(Patient{I}.MRI.post);
                    Patient{I}.MRI.reg_post = fullfile(Patient{I}.MRI.out, 'BrainPost.nii');
                    movefile(fullfile(fPath, ['r' fBase, fExt]), Patient{I}.MRI.reg_post);
                end
                if isfield(Patient{I},'CT') && isfield(Patient{I}.CT,'post')
                    [fPath, fBase, fExt] = fileparts(Patient{I}.CT.post);
                    Patient{I}.MRI.reg_postct = fullfile(Patient{I}.MRI.out, 'BrainPostCT.nii');
                    movefile(fullfile(fPath, ['r' fBase, fExt]), Patient{I}.MRI.reg_postct);
                end
                if isfield(Patient{I}.MRI,'postop')
                    [fPath, fBase, fExt] = fileparts(Patient{I}.MRI.postop);
                    Patient{I}.MRI.reg_postop = fullfile(Patient{I}.MRI.out, 'BrainPostOp.nii');
                    movefile(fullfile(fPath, ['r' fBase, fExt]), Patient{I}.MRI.reg_postop);
                end
            otherwise
                error('Not supported yet.');
        end
    end
    
        
    %% ===== CORTEX SURFACES =====
    if isfield(Patient{I}.MRI, 'pre_cortex') && ~isempty(Patient{I}.MRI.pre_cortex)
        % Canonical surfaces
        if strcmpi(Patient{I}.MRI.pre_cortex, 'canonical')
            % Computation of the SPM canonical mesh
            mesh = ImaGIN_spm_eeg_inv_mesh(Patient{I}.MRI.reg_pre, 4);
            Patient{I}.MRI.pre_cortex = mesh.tess_ctx;
        % BrainVISA folder
        elseif isdir(Patient{I}.MRI.pre_cortex)
            Patient{I}.MRI.pre_cortex = ImaGIN_load_brainvisa(Patient{I}.MRI.pre_cortex, Patient{I}.MRI.out);
        % Input surface: Single file
        elseif ischar(Patient{I}.MRI.pre_cortex)
            % Nothing specific to do: we will use directly this surface
        end
    end
    
    
    %% ===== MRI-CT FUSION =====
    % Normalize values in the CT volume, to match the amplitude of the pre-op MRI
    if strcmp(Patient{I}.MRI.ref,'pre') && isfield(Patient{I},'CT')
        % Output file name
        if isNormalize
            Patient{I}.MRI.reg_prect = fullfile(Patient{I}.MRI.out, 'wBrainPreCT.nii');
        else
            Patient{I}.MRI.reg_prect = fullfile(Patient{I}.MRI.out, 'BrainPreCT.nii');
        end
        % Load PRE
        V1 = spm_vol(Patient{I}.MRI.reg_pre);
        I1 = spm_read_vols(V1);
        % Load POST CT
        V2 = spm_vol(Patient{I}.MRI.reg_postct);
        I2 = spm_read_vols(V2);
        % Merge the two volumes
        I2 = I2./max(I2(:));
        Mask = find(I2>0.3);
        I3 = I1;
        I3(Mask) = 2*max(I1(:));%*I2(Mask);
        % Save new volume
        V3 = V1;
        V3.fname = Patient{I}.MRI.reg_prect;
        V3.dt = [8 0];
        spm_write_vol(V3,I3);
    end
end
end



%% ==================================================================================
%  ===== HELPER FUNCTIONS ===========================================================
%  ==================================================================================

%% ===== ADD COREGISTRATION =====
function matlabbatch = AddCoreg(matlabbatch, ref, source, isNormalize, outdir)
    % Add batch entry
    iReg = length(matlabbatch) + 1;
    % If normalizing: do not write output
    if isNormalize
        matlabbatch{iReg}.spm.spatial.coreg.estimate.ref      = {ref};
        matlabbatch{iReg}.spm.spatial.coreg.estimate.source   = {source};
        matlabbatch{iReg}.spm.spatial.coreg.estimate.other    = {''};
        matlabbatch{iReg}.spm.spatial.coreg.estimate.eoptions = spm_get_defaults('coreg.estimate');
    % Otherwise: write output
    else
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.ref      = {ref};
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.source   = {source};
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.other    = {''};
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.eoptions = spm_get_defaults('coreg.estimate');
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.woptions = spm_get_defaults('coreg.write');
        matlabbatch{iReg}.spm.spatial.coreg.estwrite.woptions.outdir = outdir;
    end
end 

%% ===== ADD SEGMENTATION =====
function matlabbatch = AddSegment(matlabbatch, vols)
    iSeg = length(matlabbatch) + 1;
    matlabbatch{iSeg}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{iSeg}.spm.spatial.preproc.channel.vols = {vols};
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


