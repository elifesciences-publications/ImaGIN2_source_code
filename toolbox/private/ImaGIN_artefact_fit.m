function [sig_corr,rc] = ImaGIN_artefact_fit(Signal, mmArt, flagRC, IntCor)
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
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David


%correct assuming the same RC on the series of pulses and select the RC
%value with the largest artefact power
signal=Signal(:,1:IntCor);

Uns=ones(1,size(mmArt,3));
% Power=zeros(size(signal,1),1);
if strcmp(flagRC,'same')
    R2=zeros(size(signal,1),size(mmArt,2));
     Power=inf*ones(size(signal,1),size(mmArt,2));
%     Power=zeros(size(signal,1),size(mmArt,2));
    Artefact=zeros(size(Signal,1),size(Signal,2),size(mmArt,2));
    for i1=1:size(mmArt,2)  %loop on RC
        Art=squeeze(mmArt(:,i1,:));
        m=max(abs(Art),[],2);
        index=find(m>0);
        Art=Art(index,:)./(m(index)*Uns);
        Art=unique(Art,'rows');
        
        for i2=1:size(Art,1)
            X=[Uns(1:IntCor);Art(i2,1:IntCor)];
            XX=X*X';
            
            if rcond(XX)>10^(-15)
                P=XX\X;
                beta=P*signal';
                artefactCut=beta(2,:)'*Art(i2,1:IntCor);
                Y=repmat(mean(signal), size(signal,1),1);
%                 r2=1-(sum((signal-artefactCut).^2,2))./(sum((signal-Y).^2,2));
                pow=sum(abs(artefactCut-signal),2);
                index=find(pow<Power(:,i1));
                artefact=beta(2,:)'*Art(i2,:);
                
                if ~isempty(index)
                     Power(index,i1)=pow(index);
                    Artefact(index,:,i1)=artefact(index,:);
                end
            end
            
        end
        
    end
    
    %select the RC with the best fit on all pulses
Power=mean(Power);
rc=min(find(Power==min(Power)));
    ArtefactFit=squeeze(Artefact(:,:,rc));
% [~,rc]=max(median(R2,1));
% ArtefactFit(:,:)=squeeze(Artefact(:,:,rc));
    
    
else
    
    %no constraint of same RC
    
    Power=inf*ones(size(signal,1),1);
    Artefact=zeros(size(signal),size(Signal,2));
    Art=reshape(mmArt,size(mmArt,1)*size(mmArt,2),size(mmArt,3));
    m=max(abs(Art),[],2);
    index=find(m>0);
    Art=Art(index,:)./(m(index)*Uns);
    Art=unique(Art,'rows');
    
      
    for i2=1:size(Art,1)
        
        X=[Uns(1:IntCor);Art(i2,1:IntCor)];
        P=inv(X*X')*X;
        beta=P*signal';
        artefactCut=beta(2,:)'*Art(i2,1:IntCor);
        pow=sum(abs(artefactCut-signal),2);
        index=find(pow<Power);
        artefact=beta(2,:)'*Art(i2,:);
%         pow=sum(abs(artefact),2);
%         index=find(pow>Power);
        if ~isempty(index)
            Power(index)=pow(index);
            Artefact(index,:)=artefact(index,:);
        end
    end
    
    ArtefactFit=Artefact;
    
    rc=0;
    
end


sig_corr=Signal-ArtefactFit;

end