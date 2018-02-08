function [sdy] = ImaGIN_SpikesDetection(S)
D = spm_eeg_load(S.dataset); 
t = time(D);
idx = find(t<-1);
if isempty(idx)
    idx = find(t);
end
t1 = t(idx);
x = D(:,(idx));
y = x;
%%

% First derivation estimation

%    ... bandstop filtering in 49-51 Hz band (notch filter)

f_bs_L = 49;                       % Low frequency cutoff

f_bs_H = 51;                       % High frequency cutoff

f_s     = D.fsample; 


[B,A] = butter(3,2*[f_bs_L f_bs_H]/f_s,'stop');

for n = 1:size(x,1)

    x(n,:) = filtfilt(B,A,x(n,:));

end

%% upwind third order
dy = zeros(row,col);
dy(:,1) = x(:,2)-x(:,1);
dy(:,2) = x(:,3)-x(:,2);
dy(:,3:col-2) = (2*x(:,4:col-1)+x(:,5:col)-2*x(:,2:col-3)-x(:,1:col-4))./8;
dy(:,col-1) = x(:,col) - x(:,col-1);
dy(:,col) = x(:,col) - x(:,col-1);

% Squaring
dy = (dy.^2).*sign(dy);
% Smoothing
sdy = movmean(dy,[8 8],2); 

%%

Size = 8;  % Number of channels per screenshot

if  size(D,1) < Size
    Size = size(D,1); 
    tmp = 1;
else
    tmp = floor(size(D,1)/Size);
end
for i2 = 1:tmp
    figure(i2);
    set(gcf,'Position',[629 -17 702 1101])
    for i3 = 1:Size
        subplot(Size,1,i3)
        plot(t1,(x(i3+(i2-1)*Size,:)  - mean(x(i3+(i2-1)*Size,:))  )/ max(x(i3+(i2-1)*Size,:)),  '-k','linewidth',1);hold on;
        plot(t1,(sdy(i3+(i2-1)*Size,:)- mean(sdy(i3+(i2-1)*Size,:)))/ max(sdy(i3+(i2-1)*Size,:)),'-r','linewidth',1);
        ylabel([num2str(i3+(i2-1)*Size)]);
        axis tight
    end
    zoom on
    pause;
end
    







