function [PredictionX,PredictionY] = ImaGIN_NonlinearInterdependence(x,y,TimeDelay,EmbeddingDimension,Method,TimeHorizon)
%
% REFERENCES:
%    - Schiff Phys Rev E 54(6) 6708 1997, Breakspear PhD Thesis, Quian Quiroga PhysRevE 65 041903 2002

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

if size(x,1)~=1
    x=reshape(x,1,length(x));
end
if size(y,1)~=1
    y=reshape(y,1,length(y));
end

%Normalisation
x=(x-mean(x))./std(x);
y=(y-mean(y))./std(y);

%time-delay embedding
NSample = size(x,2)-TimeDelay*(EmbeddingDimension-1)-1;
EmbeddingMatrix=[1:NSample]'*ones(1,EmbeddingDimension)+ones(NSample,1)*[0:TimeDelay:((EmbeddingDimension-1)*TimeDelay)];
X=x(EmbeddingMatrix);
Y=y(EmbeddingMatrix);

% Note: seems to be incompatible with the Theiler correction
% if exist('NSampleMax')
%     if NSample>NSampleMax
%         tmp=randperm(NSample);
%         NSample=NSampleMax;
%         X=X(tmp(1:NSample),:);
%         Y=Y(tmp(1:NSample),:);
%     end
% end
    
switch Method
    case 'Schiff'
        MeanX = mean(X);
        MeanY = mean(Y);
        %neighbourhood
        % Distance = zeros(NSample,2*EmbeddingDimension);
        NeighbourX = zeros(NSample,2*EmbeddingDimension);
        NeighbourY = zeros(NSample,2*EmbeddingDimension);
        for i1=1:NSample

            DistX=sum((ones(NSample-TimeHorizon,1)*X(i1,:)-X(1:end-TimeHorizon,:)).^2,2);
            DistY=sum((ones(NSample-TimeHorizon,1)*Y(i1,:)-Y(1:end-TimeHorizon,:)).^2,2);
            %Theiler correction
            tmp=max([1 i1-TimeDelay]):min([NSample-TimeHorizon i1+TimeDelay]);
            DistX(tmp)=0;
            DistY(tmp)=0;
            tmp=length(tmp);
            [tmp1,tmp2]=sort(DistX);
            [tmp3,tmp4]=sort(DistY);
            NeighbourX(i1,:) = tmp2(tmp+[1:2*EmbeddingDimension])';
            NeighbourY(i1,:) = tmp4(tmp+[1:2*EmbeddingDimension])';

            % [tmp1,tmp2]=sort(sum((ones(NSample,1)*X(i1,:)-X).^2,2));
            % %     Distance(i1,:) = sqrt(tmp1(2:(2*EmbeddingDimension+1)));
            % NeighbourX(i1,:) = tmp2(2:(2*EmbeddingDimension+1))';
            % [tmp1,tmp2]=sort(sum((ones(NSample,1)*Y(i1,:)-Y).^2,2));
            % %     Distance(i1,:) = sqrt(tmp1(2:(2*EmbeddingDimension+1)));
            % NeighbourY(i1,:) = tmp2(2:(2*EmbeddingDimension+1))';
        end
        %Prediction error
        PredictionXNum = zeros(NSample-TimeHorizon,1);
        PredictionXDen = zeros(NSample-TimeHorizon,1);
        PredictionYNum = zeros(NSample-TimeHorizon,1);
        PredictionYDen = zeros(NSample-TimeHorizon,1);
        for i1=1:NSample-TimeHorizon
            tmp=TimeHorizon+NeighbourX(i1,:);
            % tmp=tmp(find(tmp<=NSample));
            PredictionYNum(i1)=sqrt(sum((Y(i1+TimeHorizon,:)-mean(Y(tmp,:))).^2));
            PredictionYDen(i1)=sqrt(sum((Y(i1+TimeHorizon,:)-MeanY).^2));
            tmp=TimeHorizon+NeighbourY(i1,:);
            % tmp=tmp(find(tmp<=NSample));
            PredictionXNum(i1)=sqrt(sum((X(i1+TimeHorizon,:)-mean(X(tmp,:))).^2));
            PredictionXDen(i1)=sqrt(sum((X(i1+TimeHorizon,:)-MeanX).^2));
        end
        PredictionX=1-mean(PredictionXNum)/mean(PredictionXDen);
        PredictionY=1-mean(PredictionYNum)/mean(PredictionYDen);
        
    % case 'Quian'   %with no time horizon
    %     RXTot=zeros(1,NSample);
    %     RYTot=zeros(1,NSample);
    %     RXY=zeros(1,NSample);
    %     RYX=zeros(1,NSample);
    %     for i1=1:NSample
    %         DistX=sum((ones(NSample,1)*X(i1,:)-X).^2,2);
    %         DistY=sum((ones(NSample,1)*Y(i1,:)-Y).^2,2);
    %         %Theiler correction
    %         tmp=max([1 i1-TimeDelay]):min([NSample i1+TimeDelay]);
    %         DistX(tmp)=0;
    %         DistY(tmp)=0;
    %         tmp=length(tmp);
    %         [tmp1,tmp2]=sort(DistX);
    %         [tmp3,tmp4]=sort(DistY);
    %         RXTot(i1)=mean(tmp1(tmp+1:end));
    %         RYTot(i1)=mean(tmp3(tmp+1:end));
    %         RXY(i1)=mean(DistX(tmp4(tmp+1:(2*EmbeddingDimension+tmp))));
    %         RYX(i1)=mean(DistY(tmp2(tmp+1:(2*EmbeddingDimension+tmp))));
    %     end
    %     %Prediction error
    %     PredictionX = mean((RXTot-RXY)./RXTot);
    %     PredictionY = mean((RYTot-RYX)./RYTot);
    
    case 'Quian'   %with time horizon
        Criterion=-Inf;
        Coarse=max([floor((TimeHorizon+1)/3) 1]);
        for i0=0:Coarse:TimeHorizon
            RXTot=zeros(1,NSample-i0);
            RYTot=zeros(1,NSample-i0);
            RXY=zeros(1,NSample-i0);
            RYX=zeros(1,NSample-i0);
            for i1=1:NSample-i0
                DistX=sum((ones(NSample-i0,1)*X(i1,:)-X(1:end-i0,:)).^2,2);
                DistY=sum((ones(NSample-i0,1)*Y(i1,:)-Y(1:end-i0,:)).^2,2);
                DistXd=sum((ones(NSample-i0,1)*X(i1+i0,:)-X(1:end-i0,:)).^2,2);
                DistYd=sum((ones(NSample-i0,1)*Y(i1+i0,:)-Y(1:end-i0,:)).^2,2);
                %Theiler correction
                tmp=max([1 i1-TimeDelay]):min([NSample-i0 i1+TimeDelay]);
                DistX(tmp)=0;
                DistY(tmp)=0;
                tmp=length(tmp);
                [tmp1,tmp2]=sort(DistX);
                [tmp3,tmp4]=sort(DistY);
                RXTot(i1)=mean(tmp1(tmp+1:end));
                RYTot(i1)=mean(tmp3(tmp+1:end));
                tmptmp=i0+tmp4(tmp+1:(2*EmbeddingDimension+tmp));
                tmptmp=tmptmp(find(tmptmp<=length(DistXd)));
                RXY(i1)=mean(DistXd(tmptmp));
                tmptmp=i0+tmp2(tmp+1:(2*EmbeddingDimension+tmp));
                tmptmp=tmptmp(find(tmptmp<=length(DistXd)));
                RYX(i1)=mean(DistYd(tmptmp));
            end
            %Prediction error
            PX=mean((RXTot-RXY)./RXTot);
            PY=mean((RYTot-RYX)./RYTot);
            if Criterion<abs(PX-PY)
                Criterion=abs(PX-PY);
                PredictionX = PX;
                PredictionY = PY;
            end
        end
    
        % NeighbourX = zeros(NSample,2*EmbeddingDimension);
        % NeighbourY = zeros(NSample,2*EmbeddingDimension);
        % RYTot=zeros(1,NSample);
        % for i1=1:NSample
        %     [tmp1,tmp2]=sort(sum((ones(NSample,1)*X(i1,:)-X).^2,2));
        %     NeighbourX(i1,:) = tmp2(2:(2*EmbeddingDimension+1))';
        %     [tmp3,tmp4]=sort(sum((ones(NSample,1)*Y(i1,:)-Y).^2,2));
        %     NeighbourY(i1,:) = tmp4(2:(2*EmbeddingDimension+1))';
        %     RYTot(i1)=mean(tmp3(2:end));
        % end
        % %Prediction error
        % %         RY=zeros(1,NSample);
        % RYX=zeros(1,NSample);
        % for i1=1:NSample
        %     %             RY(i1)=mean(sum((ones(2*EmbeddingDimension,1)*Y(i1,:)-Y(NeighbourY(i1,:),:)).^2,2));
        %     RYX(i1)=mean(sum((ones(2*EmbeddingDimension,1)*Y(i1,:)-Y(NeighbourX(i1,:),:)).^2,2));
        % end
        % Prediction = mean((RYTot-RYX)./RYTot);

    case 'Chavez'
        EmbeddingMatrix2=[EmbeddingMatrix EmbeddingMatrix(:,end)+TimeHorizon];
        [C1,H1]=CorrelationIntegral([x(EmbeddingMatrix2) Y],TimeDelay,2);
        [C2,H2]=CorrelationIntegral([X Y],TimeDelay,2);
        [C3,H3]=CorrelationIntegral(x(EmbeddingMatrix2),TimeDelay,2);
        [C4,H4]=CorrelationIntegral(X,TimeDelay,2);
        PredictionX=log(C1)-log(C2)-log(C3)+log(C4);
        % PredictionX=-H1+H2+H3-H4;
        [C1,H1]=CorrelationIntegral([y(EmbeddingMatrix2) X],TimeDelay,2);
        [C2,H2]=CorrelationIntegral([X Y],TimeDelay,2);
        [C3,H3]=CorrelationIntegral(y(EmbeddingMatrix2),TimeDelay,2);
        [C4,H4]=CorrelationIntegral(Y,TimeDelay,2);
        PredictionY=log(C1)-log(C2)-log(C3)+log(C4);
        % PredictionY=-H1+H2+H3-H4;
    
end


