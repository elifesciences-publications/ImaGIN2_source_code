function mArt=ImaGIN_GenerateArtRange(rcmin, rcmax, nb, fe, art_duration, break_duration, window, mode)

id='MATLAB:NonIntegerInput';
warning('off',id);

prange=linspace(rcmin,rcmax,nb);

mArt=zeros(round(window*fe),nb);
for ii=1:length(prange)
    tau=prange(ii);
%     mArt=[mArt,ImaGIN_GenerateArtefact(tau,fe,art_duration, window, mode)];
    mArt(:,ii)=ImaGIN_GenerateArtefact(tau,fe,art_duration, break_duration, window, mode);
end

end