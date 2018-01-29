function prepare_ImaGIN_spm_eeg_filter(FilterBand, FilterType, FilterOrder, CutFreq, FileIn, FileOut)
 
% FileIn: path linking to the MEEG file to correct
% DirOut: path of the output directory
% FilterBand: band of frequency to keep ('low' or 'high')
% FilterType: type of the filter (ex: 'butterworth')
% FilterOrder: order of the filter (from 1 to...)
% CutFreq: cutoff frequency (Hz)

S = [];
S.D=FileIn;
S.filter.band = FilterBand;
S.filter.type = FilterType;

if ischar(FilterOrder)
    FilterOrder=str2num(FilterOrder);
end
S.filter.order = FilterOrder;

S.filter.dir = 'twopass';

if strcmp(FilterBand,'stop')
    CutFreq = [48 52];
end

if ischar(CutFreq)
    CutFreq=str2num(CutFreq);
end
S.filter.PHz = CutFreq;


    


%S.pre = 'f';
S.FileOut=FileOut;
D = ImaGIN_spm_eeg_filter(S);


set_final_status('OK')
end
