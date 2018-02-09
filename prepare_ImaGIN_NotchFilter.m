function prepare_ImaGIN_NotchFilter(FileIn, FileOut)

P.fname = FileIn;
P.FileOut = FileOut;
ImaGIN_NotchFilter(P);

close all

end