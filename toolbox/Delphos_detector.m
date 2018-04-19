function [results] = Delphos_detector(signal,labels, chan_type, Fs, detection_type, freq_band, artefact, thr_type,param_thr)

%%%%%

% Input:
%       signal:             nxm matrix or 1xm vector with n the number of
%                           channels and m the number of sample
%       labels:             nx1 cell array(s) of labels of your channels
%       chan_type:          string giving the type of the channel
%       Fs:                 sampling frequency
%       detection_type:     cell array(s) of string describing the type of
%                           detection; ex: {'Spk', 'Osc'}
%       freq_band:          a 1x2 matrix giving the frequency band of
%                           interested; only needed if 'Osc' selected
%       artefact:           cell arrays giving the information about the
%                           artefacted epfor i:ochs (use import_mrk for correct
%                           structure)
%       thr_type:           string, 'auto' uses the fast lfdr to find a
%                           threshold, otherwise double to set the
%                           threshold manually
%       param_thr:          Parameter of detection on the form
%                           [HFO_time_thr =  param_thr(1), HFO_freq_thr = param_thr(2).*(NbVoi/12),spike_time_thr = param_thr(3),spike_freq_thr = param_thr(4).*(NbVoi/12)]
%
%
% Ourput:
%       results:            structure => struct('markers',markers
%                           (AnyWave format),'freq_band',freq_band,
%                           'n_Spk',#Spikes,'n_Osc',#Oscillations,'labels',
%                           {labels},'detection_charac',characteristics of
%                           the events,'cfg',Input settings);
%
%
%%% Examples:
%
%   Spike detection only:       [results] = Delphos_detector(signal,labels, 'SEEG', Fs, {'Spk'}, [], [], 40, []);
%   Spike and HFO detection:    [results] = Delphos_detector(signal,labels, 'SEEG', Fs, {'Spk', 'Osc'}, [80 512], [], 40, []);
%
%%%%%


warning('off','all')
results = [];

%
% try
%     p = parcluster;
%
%     if p.NumWorkers > 6
%         parpool('local',p.NumWorkers-2);
%     else
%         parpool('local',2);
%     end
%
% end

osc_chkb_bln = any(strcmp(detection_type,'Osc')) || any(strcmp(detection_type,'Oscillation')) ...
    || any(strcmp(detection_type,'osc')) || any(strcmp(detection_type,'oscillation'));

spk_chkb_bln = any(strcmp(detection_type,'Spk')) || any(strcmp(detection_type,'Spike')) ...
    || any(strcmp(detection_type,'spk')) || any(strcmp(detection_type,'spike'));

duration = size(signal,2)/Fs;
n_sample = size(signal,2);



if isempty(thr_type)
    thr_type = 40;
elseif ischar(thr_type) && strcmp(thr_type,'auto')
    thr_type = 'auto';
end

if size(signal,1) ~= length(labels)
    errordlg('Not the same number of channels in signal and labels');
else
    n_chan = size(signal,1);
end

if ~osc_chkb_bln && ~spk_chkb_bln
    errordlg('Wrong detection type');
elseif osc_chkb_bln && size(freq_band,2) ~= 2
    errordlg('Wrong frequency band input');
end



try
    if osc_chkb_bln && ~spk_chkb_bln
        Oct(1) = floor(log(Fs./(4*max(freq_band)))/log(2));
        Oct(2) = ceil(log(Fs./(4*min(freq_band)))/log(2));
    elseif spk_chkb_bln && ~osc_chkb_bln
        Oct(1) = floor(log(Fs./(4*80))/log(2));
        Oct(2) = ceil(log(Fs./(4*8))/log(2));
    else
        Oct(1) = floor(log(Fs./(4*max(freq_band)))/log(2));
        Oct(2) = ceil(log(Fs./(4*8))/log(2));
    end
    
    if Oct(1) < -1
        Oct(1) = 0;
    end
catch ME
    errordlg(ME.message)
end

% read artefact markers
try
    if ~isempty(artefact)
        artefact_position = artefact{1,3};
        artefact_duration = artefact{1,4};
        artefact_pos = [artefact_position(:) artefact_position(:)+artefact_duration(:)];
        artefact_sample = [floor(artefact_pos(:,1)*Fs) ceil(artefact_pos(:,2)*Fs)];
        
        artefact_pos_bln = artefact_pos(:,1) <= duration & artefact_pos(:,2) <= duration;
        artefact_pos(~artefact_pos_bln,:) = [];
        
        artefact_pos_bln = artefact_pos(:,2) > duration;
        artefact_pos(artefact_pos_bln,:) = duration;
        artefact_bln = false(1,n_sample);
        
        for k = 1:size(artefact_pos,1)
            artefact_bln(artefact_sample(k,1):artefact_sample(k,2)) = true;
        end
        
    else
        artefact_pos = [];
        artefact_sample = [];
        artefact_bln = false(1,n_sample);
        
    end
catch ME
    errordlg(ME.message)
end


NbOct = length(Oct(1):Oct(2))-1;
NbVoi = 12;
VanMom = 20;

if size(param_thr,1) ~= 1 || size(param_thr,2) ~= 4
    HFO_time_thr =  1.4;
    HFO_freq_thr = 10.*(NbVoi/12);
    spike_time_thr = 1.3;
    spike_freq_thr = 11.*(NbVoi/12);
else
    HFO_time_thr =  param_thr(1);
    HFO_freq_thr = param_thr(2).*(NbVoi/12);
    spike_time_thr = param_thr(3);
    spike_freq_thr = param_thr(4).*(NbVoi/12);
end

MAX_event = cell(n_chan,1);
markers = cell(n_chan,1);

if ~strcmp(chan_type,'MEG')
    %     alpha = 10^-8;
    alpha = 0.3;
else
    alpha = 0.005;
end

[~, f] = DoG(1,Oct,NbVoi,VanMom,2,Fs,0);
n_Spk = zeros(n_chan,1);
n_Osc = zeros(n_chan,size(freq_band,1));

parfor i = 1:n_chan
    
    % As soon as a large variable, i.e tf or tfz, is not needed
    % anymore, it is set to [] to free space
    
    warning('off','all')
    %         disp ('start dog')
    % Apply DoG Analytic Continuous Wavelet Transform
    %         disp(i);
    %         disp(sum(isinf(signal(i,:))))
    %         disp(sum(isnan(signal(i,:))))
    [tf, ~] = DoG(signal(i,:),Oct,NbVoi,VanMom,2,Fs,0);
    %         disp(sum(isnan(tf(:))))
    %         disp(sum(isinf(tf(:))))
    %         disp(isempty(tf))
    %         disp(sum(isinf(tf(:))))
    %         disp ('hello')
    % Normalized the TF
    [tfz, ~, tf_real_modif, sigma] = z_H0( tf, Fs, artefact_bln);
    tf = [];
    
    if strcmp(thr_type, 'auto')
        
        tf_stat_HFO = tf_real_modif(~artefact_bln,:);
        
        N = length(tf_stat_HFO);
        
        if N > 16000;
            decimate = floor(linspace(Fs,N-Fs,15000));
        else
            decimate = Fs:N-Fs;
        end
        if (NbOct*NbVoi - 3*NbVoi) > 0 && osc_chkb_bln && f(end) == Fs/4 && Fs > 1000;
            
            tf_stat_low = tf_stat_HFO(decimate,1:(NbOct*NbVoi - 3*NbVoi-1));
            tf_stat_HFO = tf_stat_HFO(decimate,(NbOct*NbVoi - 3*NbVoi):(NbOct*NbVoi));
            [thr_low1, thr_high1] = fast_lfdr(tf_stat_HFO(:), alpha);
            [thr_low2, thr_high2] = fast_lfdr(tf_stat_low(:), alpha);
            tf_stat_HFO = [];
            tf_stat_low = [];
            
            thr_low = max([thr_low1, thr_low2]); %max since negative
            thr_high = min([thr_high1, thr_high2]);
        else
            tf_stat_HFO = tf_stat_HFO(decimate,:);
            [thr_low, thr_high] = fast_lfdr(tf_stat_HFO(:), alpha);
            tf_stat_HFO = [];
        end
        %                 tf_stat = tf_real_modif(decimate,(NbOct*NbVoi - 3*NbVoi):(NbOct*NbVoi));
        tf_real_modif = [];
        %             [thr_low, thr_high] = lfdr(tf_stat(:), 0.03, 'tLocationScale', 0);
        
        %                 [thr_low, thr_high] = fast_lfdr(tf_stat(:), alpha);
        thr = single(mean([-thr_low thr_high])^2);
        tf_z_thr = tfz.*(tfz > thr);
        %                 tf_stat = [];
    else
        tf_real_modif = [];
        thr = single(thr_type);
        tf_z_thr = tfz.*(tfz > thr);
    end
    [ t_max,f_max,value_max ] = maxima_tf(tf_z_thr);
    tf_z_thr = [];
    
    if ~isempty(t_max)
        
        [t_max, i_sort] = sort(t_max);
        f_max = f_max(i_sort);
        value_max = value_max(i_sort);
        [ DL, DH, Area, L1, L2, H1, H2] = region_half_high_charac( tfz, t_max, f_max, value_max, Fs );
        tfz = [];
        
        FWHMt = ceil(1/pi*sqrt(2*log(2)*VanMom)*Fs./f(f_max));
        % FWHMf = ceil(2*pi*sqrt(log(2)/(2*VanMom))*f(f_max)');
        
        MAX = [t_max, f_max, value_max, DL, DH, Area, L1, L2, H1, H2, FWHMt];
        
        
        freq_thr = find(f > 120, 1, 'first');
        if isempty(freq_thr)
            freq_thr = max(MAX(:,2)) + 1;
        end
        
        
        MAX_HFO = [];
        MAX_IES = [];
        
        if spk_chkb_bln
            
            spike = MAX(:,4)./MAX(:,end) < spike_time_thr & MAX(:,5) >= spike_freq_thr & MAX(:,2) < freq_thr...
                & ~artefact_bln(MAX(:,1))';
            
            [ MAX_IES, ~] = sparse_detection_selection( MAX(spike,:),Fs,sigma,'Spk');
        end
        
        if osc_chkb_bln
            oscill = MAX(:,4)./MAX(:,end) >= HFO_time_thr & MAX(:,5) < HFO_freq_thr & ...
                f(MAX(:,2)) >= min(freq_band) & f(MAX(:,2)) <= max(freq_band) & ...
                ~artefact_bln(MAX(:,1))';
            [ MAX_HFO, ~] = sparse_detection_selection( MAX(oscill,:),Fs,sigma,'Osc');
        end
        
        
        sigma = [];
        
        n_Osc(i) = size(MAX_HFO,1);
        n_Spk(i) = size(MAX_IES,1);
        MAX_event{i} = [MAX_HFO; MAX_IES]';
        
        
        
        %             if osc_chkb_bln
        %                 freq_oscill = f(MAX_HFO(:,2)); % extract the frequency of the oscillations
        %                 for k = 1:size(freq_band,1)
        %
        %                     bln = (freq_oscill >= freq_band(k,1)) & (freq_oscill < freq_band(k,2)); % test if it is in this freq band
        %
        %                     n_Osc(i,k) = sum(bln); % compute rate, duration is in sec
        %
        %                 end
        %             end
        
        temp_handles = struct('n_Osc',size(MAX_HFO,1),'n_Spk',n_Spk(i),'MAX_event',[MAX_HFO; MAX_IES],...
            'Fs',Fs,'f',f,'selected_channel',i,'labels',{labels});
        [markers{i}] = awt_detection2marker(temp_handles);
        
    end
    
    %         temp_handles = struct('markers',markers{i},'freq_band',freq_band,'rate_spk',...
    %             rate_spk,'rate_osc',rate_osc,'spk_display',spk_display,'labels',{labels},'duration',duration);
    %
    
end

MAX_event = [MAX_event{:}]';

markers = [markers{:}];

cfg = struct('Fs',Fs,'f',f,'NbVoi',NbVoi,'NbOct',NbOct,...
    'VanMom',VanMom,'duration',duration);
results = struct('markers',markers,'freq_band',freq_band,'n_Spk',n_Spk,'n_Osc',n_Osc,'labels',{labels},'detection_charac',MAX_event,'cfg',cfg);

% catch ME
%
%     errordlg(ME.message)
%
% end

function [wt, freqlist] = DoG(sig,Oct,NbVoi,VanMom,Nexp,Fs,scaling)
% CWT:	Continuous wavelet transform
% usage:	wt = cwt(sig,Oct,NbVoi,VanMom)

sig = sig(:);
siglength = length(sig);
fff = (0:(siglength-1))*2*pi/siglength;

Cst = 4*VanMom/(pi*pi);
fsig = fft(sig);
% disp(sum(isinf(fsig(:))))
% disp(sum(isnan(fsig(:))))
NbOct = length(Oct(1):Oct(2))-1;

wt = complex(zeros(NbOct*NbVoi,siglength,'single'));
freqlist=zeros(NbOct*NbVoi,1,'single');

j=1;
for oct = Oct(1):(Oct(2)-1)
    for voi = 0:(NbVoi-1)
        scale = 2^(oct + voi/NbVoi);
        freqlist(j)=Fs/(4*scale);
        %         tmp = double(scale * fff);
        %         psi = double((tmp.^VanMom)) .* exp(-Cst*tmp.^Nexp/2);
        tmp = scale * fff;
        psi = (tmp.^VanMom).* exp(-Cst*tmp.^Nexp/2);
        fTrans = fsig .* psi';
        if scaling
            wt(j,:) = ifft(fTrans)';
        else
            wt(j,:) = sqrt(scale)*ifft(fTrans)';
        end
        j=j+1;
    end
end

wt = flipud(wt)';
freqlist = flip(freqlist);




function [thr_low, thr_high] = fast_lfdr(y, alpha)

b = figure('Visible','off');
h = histogram(y,'Normalization','pdf');
f = h.Values;
x = mean([h.BinEdges(2:end); h.BinEdges(1:end-1)]);
f0 = normpdf(x,0,1);
lf = log(f0./f);
thr_low = x(find(lf<log(alpha) & x < 0,1,'last'));
thr_high = x(find(lf<log(alpha) & x > 0,1,'first'));
if isempty(thr_low) && isempty(thr_high)
    thr_low = -Inf;
    thr_high = Inf;
elseif isempty(thr_low) && ~isempty(thr_high)
    thr_low = -thr_high;
elseif ~isempty(thr_low) && isempty(thr_high)
    thr_high = -thr_low;
end
close(b);
function [ markers ] = awt_detection2marker( handles )

N = size(handles.MAX_event,1);


labels_freqband = {'Very Fast Osc' 'Fast Ripple' 'Ripple'...
    'Gamma' 'Beta' 'Alpha' 'Theta' 'Delta' 'Infra slow'};
freqband = [ 500,Inf; 250,500; 80,250;...
    24,80; 12.4,24; 7.4,12.4; 3.5,7.4; 1,3.5; 0,1];
freqcolor = {'#ff00ff' '#ff0000' '#ff8000' '#ffb000'...
    '#00b000' '#00b0b0' '#0070ff' '#0000ff' '#000090'};

label_spk = 'Spike';
color_spk = '#303030';
default_name = 'Oscillation';
default_color = '#c0c0c0';

% markers = struct('label','','value',0,'position',0,'duration',0,'channels',{});
%
% for i = 1:N
%
%     if i <= handles.n_Osc
%         markers(1,i).label = ['TFD Osc: ' num2str(handles.f(handles.MAX_event(i,2)),3) 'Hz'];
%         markers(1,i).value = (handles.f(handles.MAX_event(i,2)) >= 250)*2 + ...
%         (handles.f(handles.MAX_event(i,2)) >= 80 &  handles.f(handles.MAX_event(i,2)) < 250)*1;
%     else
%         markers(1,i).label = ['TFD Spk: ' num2str(handles.f(handles.MAX_event(i,2)),3) 'Hz,'...
%             ' Sup = ' num2str(handles.MAX_event(i,6),3)];
%         markers(1,i).value = -1;
%     end
%
%
%     markers(1,i).position = double(handles.MAX_event(i,1)/handles.Fs);
%     markers(1,i).duration = 0;
%     markers(1,i).channels = handles.labels(handles.selected_channel);
%
% end

if N > 0
    L = repmat({default_name},handles.n_Osc,1);
    C = repmat({default_color},handles.n_Osc,1);
    for i = 1:size(freqband,1)
        
        bln = handles.f(handles.MAX_event(1:handles.n_Osc,2)) >= freqband(i,1) & ...
            handles.f(handles.MAX_event(1:handles.n_Osc,2)) < freqband(i,2);
        
        L(bln) =  labels_freqband(i);
        C(bln) =  freqcolor(i);
        
    end
    
    L = [L; repmat({label_spk},N-handles.n_Osc,1)];
    C = [C; repmat({color_spk},N-handles.n_Osc,1)];
    
    
    markers = struct('label',L,...
        'value',num2cell(double(handles.f(handles.MAX_event(:,2)))),...
        'position',num2cell(double(handles.MAX_event(:,1)/handles.Fs)),'duration',0,...
        'channels',repmat({handles.labels(handles.selected_channel)},N,1),...
        'color', C)';
    
    %     save(fullfile(handles.infos.plugin_dir,'markers.mat'),'markers')
else
    
    markers = [];
    
end

function [ t_max,f_max,value_max ] = maxima_tf( image )
% Find the local maxima (8-connectivity)in time-frequency space and gives back its
% amplitudes and its frequency (/!\ this function is made for (t,f) space
% generated using awt_freqlist, which means the input image is the
% transpose of the normal(t,f) space)
%
% Input:
%
%   image:      time-frequency space (works for anytime of image actally)
%               in which you want to find the local maxima
%
% Outputs:
%
%   t_max:      lines of the local maxima, which refers to the time axis of
%               the (t,f) space
%   f_max:      column of the local maxima, which refers to the frequency
%               axis of the (t,f) space
%   value_max:  values of the local maxima
%
% Author: Nicolas Roehri (INS Marseille)


image(isnan(image)) = 0;

N = size(image,1);
M = size(image,2);

if M*N < 5*10^5
    
    M = 2*max(image(:));
    
    horizontal_test1 = (image - [ M*ones(size(image,1),1) image(:,1:end-1)] ) >=0;
    horizontal_test2 = (image - [ image(:,2:end) M*ones(size(image,1),1)]) >=0;
    
    vertical_test1 = (image - [ M*ones(1,size(image,2)); image(1:end-1,:)]) >=0;
    vertical_test2 = (image - [ image(2:end,:); M*ones(1,size(image,2))]) >=0;
    
    diagonal_test1 = (image - [ M*ones(size(image,1)-1,1) image(2:end,1:end-1); ...
        M*ones(1,size(image,2))]) >=0;
    diagonal_test2 = (image - [  image(2:end,2:end) M*ones(size(image,1)-1,1); ...
        M*ones(1,size(image,2))]) >=0;
    
    diagonal_test3 = (image - [  M*ones(1,size(image,2)); ...
        M*ones(size(image,1)-1,1) image(1:end-1,1:end-1)]) >=0;
    diagonal_test4 = (image - [  M*ones(1,size(image,2)); ...
        image(1:end-1,2:end) M*ones(size(image,1)-1,1)]) >=0;
    
    
    % result_test = 1/4*(horizontal_test1 + horizontal_test2 + vertical_test1 + vertical_test2) == 1;
    result_test = 1/8*(horizontal_test1 + horizontal_test2 + vertical_test1 + vertical_test2 ...
        + diagonal_test1 + diagonal_test2 + diagonal_test3 + diagonal_test4) == 1;
    
    % maximum = max(max(result_test));
    % result_test(result_test < maximum*0.5) = 0;
    
else
    
    mask1 = [-1   0   0
        0   1   0
        0   0   0];
    
    mask2 = [ 0  -1   0
        0   1   0
        0   0   0];
    
   
    result_test = true(size(image));
    
    for i = 1:4
        
        result_test = result_test.*(conv2(image,mask2,'same') >= 0);
        result_test = result_test.*(conv2(image,mask1,'same') >= 0);
        mask2 = rot90(mask2);
        mask1 = rot90(mask1);
        
    end
    
    result_test(1:end,[1,M]) = 0;
    result_test([1,N],1:end) = 0;
    
end

result_test = result_test.*image;


[t_max,f_max,value_max] = find(result_test);



function [ DL, DH, Area, L1, L2, H1, H2, tH1, tH2 ] = region_half_high_charac( tf_z_thr, t_max, f_max, value_max, Fs )
% This function measures the FWHM in time and in frequency from the local maxima
%
% The function finds the first value below the half of the local maximum or
% the one above the local maximum in case there is two local maxima colliding
%

H1 = zeros(size(t_max));
H2 = zeros(size(t_max));
L1 = zeros(size(t_max));
L2 = zeros(size(t_max));

for k = 1:length(t_max)
    
    T = 0.5*Fs;
    
    if t_max(k)+T > size(tf_z_thr,1)
        T = size(tf_z_thr,1)-t_max(k);
    elseif t_max(k) <= T
        T = t_max(k)-1;
    end
    
    a = find(tf_z_thr(t_max(k),f_max(k)+1:end) < value_max(k)/2 | tf_z_thr(t_max(k),f_max(k)+1:end) > value_max(k),1,'first'); % f^^^^
    b = f_max(k) - find(tf_z_thr(t_max(k),1:f_max(k)) < value_max(k)/2 | tf_z_thr(t_max(k),1:f_max(k)) > value_max(k),1,'last'); % f\/\/\/ 
    c = find(tf_z_thr(t_max(k)+1:t_max(k)+T,f_max(k)) < value_max(k)/2 | tf_z_thr(t_max(k)+1:t_max(k)+T,f_max(k)) > value_max(k) ,1,'first'); % t===>
    d = (T+1) - find(tf_z_thr(t_max(k)-T:t_max(k),f_max(k)) < value_max(k)/2 | tf_z_thr(t_max(k)-T:t_max(k),f_max(k)) > value_max(k),1,'last'); % <===t 
   
    % Boundary conditions
    % In frequency space, symmetry is preferred
    if isempty(a)
        if isempty(b)
            H1(k,1) = size(tf_z_thr,2)-f_max(k);
        else
            H1(k,1) = b;
        end
    else
        H1(k,1) = a;
    end
    
    if isempty(b)
        if isempty(a)
            H2(k,1) = f_max(k)-1;
        else
            H2(k,1) = a;
        end
    else
        H2(k,1) = b;
    end
    
    if isempty(c)
        L1(k,1) = T/2;   
    else
        L1(k,1) = c;
    end
    
    if isempty(d)
        L2(k,1) = T/2;
    else
        L2(k,1) = d;
    end
    
    
end
DL = L1+L2;
DH = H1+H2;
Area = pi.*DH.*DL/4;

function [ MAX, b] = sparse_detection_selection(MAX,Fs,sigma,type)
%%%% MAX_feature = [t_max, f_max, value_max, DL, DH, Area, L1, L2, H1, H2, FWHMt];
N = size(MAX,1);
b = true(N,1);

%%Two intervals i1 = (s1, e1) and i2 = (s2, e2) overlap if and only if s2 <
%%e1 and s1 < e2 http://www.rgrjr.com/emacs/overlap.html

if strcmp(type,'Osc') 
    a = 2;
elseif strcmp(type,'Spk')
    a = 1;
end


for i = 1:N
    
    ts1 = MAX(i,1)-a*MAX(i,8);
    te1 = MAX(i,1)+a*MAX(i,7);
    
    fs1 = MAX(i,2)-MAX(i,10);
    fe1 = MAX(i,2)+MAX(i,9);
    
    S = diff(abs(MAX(:,1)-MAX(i,1))) < Fs;
    
    s = find(S, 1, 'first');
    e = find(S, 1, 'last');
    
    for j = [s:i-1 i+1:e]
        

        ts2 = MAX(j,1)-a*MAX(j,8);
        te2 = MAX(j,1)+a*MAX(j,7);
        
        fs2 =  MAX(j,2)-MAX(j,10);
        fe2 = MAX(j,2)+MAX(j,9);

        
        if strcmp(type,'Osc')
            
            if (ts2<=te1 && ts1<=te2) && (fs2<=fe1 && fs1<=fe2) %% if overlap in time and freq => choose the more powerful
                
                if MAX(j,3) >= MAX(i,3)
                    
                    b(i) = false;
                    break
                end
            elseif (ts2<=te1 && ts1<=te2) && (12+MAX(j,2)-3<MAX(i,2) && MAX(i,2)<12+MAX(j,2)+3) %% if overlap in time and twice the freq => choose the lower freq (harmonic issue)
                if MAX(j,3)*(sigma(MAX(j,2))).^2 >= MAX(i,3)*(sigma(MAX(i,2))).^2
                    
                    b(i) = false;
                    break
                end
                break
            end
        elseif strcmp(type,'Spk') && (ts2<=te1 && ts1<=te2)
            
            if MAX(i,2) >= MAX(j,2)
                
                b(i) = false;
                break
            end
            
        end
        
        
    end
    
end

MAX = MAX(b,:);


function [ tf_z, tf2, tf_real_modif, sigma_real_N, Q] = z_H0( tf, Fs, artefact_bln)
if size(tf,2) > size(tf,1)
    tf = tf';
end

w = single(tukeywin(size(tf,1),0.25*Fs/size(tf,1))*ones(1,size(tf,2)));
% tf = w.*tf;

tf_real_modif = real(tf);
tf_imag_modif = imag(tf);
tf2 = abs(tf).^2;
tf = [];
Nf = size(tf_real_modif,2);
% Nt = size(tf,1);
% mu_real_N = zeros(1,Nf);
sigma_real_N = zeros(1,Nf,'single');
Q = 0;

% reject = check_poisson(tf2,Fs);

% tf_stat = tf_real(~reject,:);
tf_stat = tf_real_modif;
tf_stat(artefact_bln,:) = []; % reject part labeled as artefact
N = size(tf_stat,1);

if N > 16000;
    decimate = floor(linspace(Fs,N-Fs,15000));
else
    decimate = Fs:N-Fs;
end
tf_stat = tf_stat(decimate,:);
K = 1.5;
% f0 = zeros(Nf,40-1);
% C = cell(1,Nf);

% Mu = [0; 0];
% Sigma(:,:,1) = 1;
% Sigma(:,:,2) = 10;
% PComponents = [0.3,0.7];
% options = statset('MaxIter',150);


for i = 1:Nf
    b = tf_stat(:,i);
    
    % b = mean(repmat(b,1,500) + 2*max(b)*(rand(length(b),500)-0.5),2);
    IQR = iqr(b);
    q = quantile(b, [1/4 3/4]) + [-K*IQR K*IQR];
    b(b<q(1) | b>q(2)) = [];
        pd = fitdist(b,'Normal');
    
    tf_real_modif(:,i) = tf_real_modif(:,i)/pd.sigma;
    tf_imag_modif(:,i) = tf_imag_modif(:,i)/pd.sigma;
    sigma_real_N(1,i) = pd.sigma;
    
    % C{1,i} = pd;
end
b = [];
tf_stat = [];
tf_z = w.*abs(tf_real_modif + 1i*tf_imag_modif).^2;
tf_imag_modif = [];
% tfz_name = ['006GRE_TFZ_' chl '.mat'];
% tfz_name = strrep(tfz_name,'''','p');
% save(tfz_name,'tf_z');





