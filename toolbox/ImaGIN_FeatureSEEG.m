function T =  ImaGIN_FeatureSEEG(S)

sFile = S.FileName;
[pth, fName, ~] = fileparts(sFile);
clear S
S.Fname = sFile;
S.EventType = 'Stim';
%S.StartInterpolation= -0.015;
%S.EndInterpolation  = 0.015;
S.StartInterpolation = -0.008;
S.EndInterpolation   = 0.008;

D = ImaGIN_InterpolationFilter(S);
clear S

iFile = [pth '/i' fName];
nFile = [pth '/ni' fName];
P.Fname =  iFile;
P.Freq  = 50;
ImaGIN_NotchFilter(P) % notch filter 50Hz

P.Fname =  nFile;
P.Freq  = 100;
ImaGIN_NotchFilter(P) % notch filter 100Hz
clear P.Freq;
delete([nFile '.*']);
delete([iFile '.*']);
nnFile  =  [pth '/nni' fName];
P.Fname =  nnFile;
P.Freq  = 150;
ImaGIN_NotchFilter(P) % notch filter 150Hz
clear P.Freq;
delete([nnFile,'.*']);

nnnFile  =  [pth '/nnni' fName];
P.LFname = nnnFile;
ImaGIN_LowPassFilter(P) % lowpass filter 0.2Hz
lpf_nFile = [pth '/lpf_nnni' fName];
delete([nnnFile,'.*'])

D = spm_eeg_load(lpf_nFile);
sens= indchantype(D,'eeg');
elec= sensors(D,'eeg');
pos = elec.elecpos; 
nx  = size(sens,2);
nn  = 10;
nt  = find((time(D)>= -0.5));
ny  = D.nsamples;
data= D(1:nx,nt(1):ny);
logScale = 1;
%% 
rawVar  = var(data, [], 2);       % Compute raw data variance
ch_mean = mean(data,2);           
ch_ampl = range(data,2);
ch_grad = mean(diff(data,1,2),2); % Channel median gradient
ch_kurt = kurtosis(data,1,2);     % kurtosis

rawVar  = nonzero_noninf(rawVar); % Check that there is no zero variance channels
ch_mean = nonzero_noninf(ch_mean);
ch_ampl = nonzero_noninf(ch_ampl);
ch_grad = nonzero_noninf(ch_grad);
ch_kurt = nonzero_noninf(ch_kurt);

bandVar = rawVar;
ch_dev  = ch_mean - mean(ch_mean); % deviation 

ch_xcorr= zeros(nx, 1);   
noIdx  = zeros(nx, 1); 
ch_hurs  = zeros(nx, 1);

m_var  = zeros(nx, 1);

for i = 1:nx
    d = misc_euclidean_dist(pos, pos(i,:)); % Euclidean distance between channels
    [~, sortIdx] = sort(d, 'ascend');
    idx = sortIdx(2:min(nn, numel(d))); 
    
    ch_hurs(i)= hurst(data(i,:));           % Computer hurst exponent of time series
    m_var(i)  = median(bandVar(sortIdx(1:min(nn, numel(d))))); 
 
 
    for j = 1:numel(idx)
        thisCorr = xcorr(data(i,:), data(idx(j),:), 0, 'coeff');
        if isnan(thisCorr)
            % NaN due to one of either of the two channels being flat
            continue;
        end
        ch_xcorr(i) = ch_xcorr(i) + abs(thisCorr);
    end
    noIdx(i) = i;
    ch_xcorr(i) = ch_xcorr(i)/numel(idx); 
end

m_var(m_var   < 1e-3) = median(m_var);

% Normalize based on local values
if ~isempty(nn) && ~isinf(nn)
    bandVar = bandVar./m_var;

end
if logScale 
    ch_var = 10*log10(bandVar);
else
    ch_var = bandVar;
end
%%

T = table(noIdx, ch_xcorr, ch_var, ch_dev, ch_ampl, ch_grad, ch_kurt, ch_hurs);

% csvfilename = ['-csv_' fName '.csv'];
% writetable(T,csvfilename,'Delimiter',',');

delete([lpf_nFile,'.*']);

end


function bandVar = nonzero_noninf(bandVar)

isZero = (abs(bandVar) < eps);

if all(isZero)
    error('bad_channels:ZeroValue', ...
        'All data channels have zero-value!');
end

bandVar(isZero) = 1e-3*min(bandVar(~isZero));
isInf = (bandVar > 1e1000);
isInf(bandVar(isInf) < 1000*max(bandVar(~isInf))) = false;

if all(isInf)
    error('bad_channels:InfValue', ...
        'All data channels have infinite value!');
end
bandVar(isInf) = 1e3*max(bandVar(~isInf));
end



function h =  hurst(x)
%      Estimate Hurst exponent on a timeseries.
%  
%      The estimation is based on the second order discrete derivative.
%  
%      Parameters
%      ----------
%      x : 1D array
%          The timeseries to estimate the Hurst exponent for.
%  
%      Returns
%      -------
%      h : float
%          The estimation of the Hurst exponent for the given timeseries.
     
y = cumsum(diff(x,1));

b1 = [1, -2, 1];
b2 = [1,  0, -2, 0, 1];

% Second order derivative
y1 = filter(b1, 1, y);
% First values contain filter artifacts
y1 = y1(length(b1)+1:end-1);
% Wider second order derivative
y2 = filter(b2, 1, y);
% First values contain filter artifacts
y2 = y2(length(b2)+1:end-1);
s1 = mean(y1 .^2);
s2 = mean(y2 .^2);

h = 0.5*log2(s2/s1); % Hurst exponent index

end

function y = misc_euclidean_dist(a,b)
% EUCLIDEAN_DIST - Euclidean distance between two points
%
%   EUCLIDEANDIST(A,B) where A and B give the three cartesian coordinates
%   of two points in the space returns the Euclidean distance between those
%   two points.
%
%   EUCLIDEANDIST(A,B) where A and B are two matrices of dimensions Kx3
%   returns a vector of dimensions Kx1 with the Euclidean distance between
%   the rows of A and the rows of B. 
%
%   EUCLIDEANDIST(A,B) where A (or B) is a vector and the other input
%   is a Kx3 matrix returns the Euclidean distance between the single point
%   in the space given by A (or B) and the multiple points given by the
%   other input parameter.

if nargin < 2, 
    ME = MException('euclideanDistance:needMoreInputs','Some input parameters are missing');
    throw(ME);
end

if size(a,1) ~= size(b,1) && size(a,1)>1 && size(b,1) > 1,
    ME = MException('euclideanDist:invalidDim','Invalid dimensions in input parameter.');
    throw(ME);
end

if size(a,1) > size(b,1),
    tmp = a;
    a = b;
    b = tmp;
    clear tmp;
elseif size(a,1)>1 && size(a,1) == size(b,1),
    y = nan(size(a,1),1);
    for i = 1:size(a,1)
       y(i) = misc_euclidean_dist(a(i,:), b(i,:));       
    end
    return;
end
   
a = repmat(a,size(b,1),1);
y = sqrt(sum((a-b).^2,2));
end
