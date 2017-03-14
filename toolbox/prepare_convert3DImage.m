function prepare_convert3DImage(CorticalMesh, SaveMNI, pathMRI, FileIn, FileOut)


% CorticalMesh: 0 or 1, indicates whether a cortical mesh must be used or
% not
% sMRI: path to the normalised MRI
% SaveMNI: indicates whether to save in the native landmark or in MNI


    clear SS
    SS.SizeSphere=5;
    SS.SizeHorizon=10;
    SS.n=3;
    SS.TimeWindow=0;
    SS.TimeWindowWidth=0;
    SS.interpolate_bad=0;
    SS.CorticalMesh=CorticalMesh;
    SS.sMRI=pathMRI;
    SS.Atlas='Human';
    SS.SaveMNI=SaveMNI;
    SS.FileOut=FileOut;
    
    if isdir(FileIn)
        cd(FileIn)
        rmdir(FileIn,'s')
    end
    
    SS.Fname=FileIn;
    
    try
        ImaGIN_spm_eeg_convertmat2ana_3D(SS)
    end
    
    %smooth image
    if ~CorticalMesh
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {FileIn,'sample_0.img,1'};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 1;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm('defaults', 'EEG');
        spm_jobman('run', matlabbatch);
    end
    
end

