function prepare_MRImplantation(ElectrodesTxtFile, FileSource, DeformationFile, DirOut)

%Add electrode position
clear S
S.FileSource=FileSource;
S.FileName=ElectrodesTxtFile;
S.Deformation=DeformationFile;
S.DirOut=DirOut;
ImaGIN_MRImplantation(S);

end