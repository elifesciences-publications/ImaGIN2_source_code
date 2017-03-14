function [Cel, Cind, x, y, z, Ic] = ImaGIN_spm_eeg_locate_channels(D, n, interpolate_bad,dmax,Ctf,Atlas)
% function [Cel, Cind, x, y, z] = ImaGIN_spm_eeg_locate_channels(D, n, interpolate_bad)
%
% Locates channels and generates mask for converting EEG data to analyze
% format on the scalp
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_locate_channels.m 112 2005-05-04 18:20:52Z john $

switch Atlas
    case{'Human'}
        bb = [[-78 -112 -50];[78 76 85]];
        tmp=spm('Defaults','EEG');
        bb=tmp.normalise.write.bb;
    case{'Rat'}
        bb = [[-80 -156 -120];[80 60 10]];
    case{'Mouse'}
        bb = [[-48 -94 -70];[48 72 0]];
    case{'PPN'}
        bb = [[-8 -5 -20];[8 6 2]];     %Brainstem full
end


% find channel positions in 3D
%if size(D.channels.order,2)==2
    Nchannels=(1+sqrt(1+8*D.nchannels))/2;
    Nchannels=D.nchannels;
    Nchannels=length(Ctf.label);
%     D.channels.eeg=1:Nchannels;
%end
    
try
    if ~isfield(Ctf,'pnt')
        Ctf.pnt=Ctf.elecpos;
    end
%     Cel = Ctf.Cpos2(:, D.channels.order(D.channels.eeg));
    Cel = Ctf.pnt';    
catch
    Nchannels=(1+sqrt(1+8*D.nchannels))/2;
    Nchannels=D.nchannels;
    Nchannels=length(Ctf.label);
    D.channels.order=1:Nchannels;
    D.channels.eeg=1:Nchannels;
%     Cel = Ctf.Cpos2(:, D.channels.order(D.channels.eeg));
    Cel = Ctf.pnt';    
end

% Bad = [];
% if isfield(D, 'badchannels')
%     Bad = D.badchannels(:);
% end
Bad=badchannels(D);

% For mapping indices
Itmp = zeros(1,Nchannels);
Itmp(1:Nchannels) = 1:Nchannels;        %to be optimised ?

Cind = setdiff(1:Nchannels, Bad);

[x,y,z] = meshgrid(bb(1,1):n:bb(2,1),...
    bb(1,2):n:bb(2,2),...
    bb(1,3):n:bb(2,3));

% dmax=[];
% for i1=1:size(Cel,2)
%     i2=setdiff(1:size(Cel,2),i1);
%     tmp=(Cel(1,i2)-Cel(1,i1)).^2+(Cel(2,i2)-Cel(2,i1)).^2+(Cel(3,i2)-Cel(3,i1)).^2;
%     dmax=[dmax min(tmp)];
% end
% dmax=max(dmax);
% % dmax=mean([max(dmax) min(dmax)]);
dmax=dmax^2;  %5mm

if interpolate_bad
    % keep bad electrode positions in
    Ic=[];
    for i1=1:size(Cel,2)
        d=(x(:)-Cel(1,i1)).^2+(y(:)-Cel(2,i1)).^2+(z(:)-Cel(3,i1)).^2;
        Ic=[Ic;find(d<=dmax)];
    end
else
    % or don't
    Ic=[];
    tmp=Itmp(Cind);
    for i1=1:length(Cind)
        d=(x(:)-Cel(1,tmp(i1))).^2+(y(:)-Cel(2,tmp(i1))).^2+(z(:)-Cel(3,tmp(i1))).^2;
        Ic=[Ic;find(d<=dmax)];
    end
end
Ic=unique(Ic);

% %erosion
% ninit=ceil(5/n);
% ninit=1;
% tmp=zeros(size(x));
% tmp(Ic)=1;
% [a,b,c]=ind2sub(size(x),Ic);
% % for i1=1:ninit
% %     for i2=1:length(Ic)
% %         if ~(tmp(min([size(tmp,1) a(i2)+1]),b(i2),c(i2))==1&...
% %                 tmp(max([1 a(i2)-1]),b(i2),c(i2))==1&...
% %                 tmp(a(i2),min([size(tmp,2) b(i2)+1]),c(i2))==1&...
% %                 tmp(a(i2),max([1 b(i2)-1]),c(i2))==1&...
% %                 tmp(a(i2),b(i2),min([size(tmp,3) c(i2)+1]))==1&...
% %                 tmp(a(i2),b(i2),max([1 c(i2)-1]))==1)
% %             Ic(i2)=0;
% %         end
% %     end
% %     Ic=Ic(find(Ic));
% %     [a,b,c]=ind2sub(size(x),Ic);
% % end
% for i1=1:ninit
%     for i2=1:length(Ic)
%         if ~(tmp(min([size(tmp,1) a(i2)+1]),b(i2),c(i2))==1&...
%                 tmp(max([1 a(i2)-1]),b(i2),c(i2))==1)
%             Ic(i2)=0;
%         end
%     end
%     Ic=Ic(find(Ic));
%     [a,b,c]=ind2sub(size(x),Ic);
%     for i2=1:length(Ic)
%         if ~(tmp(a(i2),min([size(tmp,2) b(i2)+1]),c(i2))==1&...
%                 tmp(a(i2),max([1 b(i2)-1]),c(i2))==1)
%             Ic(i2)=0;
%         end
%     end
%     Ic=Ic(find(Ic));
%     [a,b,c]=ind2sub(size(x),Ic);
%     for i2=1:length(Ic)
%         if ~(tmp(a(i2),b(i2),min([size(tmp,3) c(i2)+1]))==1&...
%                 tmp(a(i2),b(i2),max([1 c(i2)-1]))==1)
%             Ic(i2)=0;
%         end
%     end
%     Ic=Ic(find(Ic));
%     [a,b,c]=ind2sub(size(x),Ic);
% end
    
            

x = x(Ic); y = y(Ic); z = z(Ic);

Cel=Cel';
Cel = Cel(Itmp(Cind), :);
if length(Cind)==D.nchannels
    Cind=1:D.nchannels;
end
