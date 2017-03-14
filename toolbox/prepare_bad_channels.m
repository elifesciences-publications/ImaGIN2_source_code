function prepare_bad_channels(Badp, FileIn, FileOut)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% Badp: path linking to the file with the name of bad channels

Bad=load(Badp);
[DirOut,NameOut,~]=fileparts(FileOut);
if ~isempty(Bad)
    S.D=FileIn;
    D=spm_eeg_load(S.D);
    D=badchannels(D, Bad, 1);
    nameDproc=fullfile(DirOut, NameOut);
    D2=clone(D,nameDproc, [D.nchannels D.nsamples D.ntrials]);
    D2(:,:,:)=D(:,:,:);
    save(D);
end

end