function [S,COH]=ImaGIN_mvfreqz_SCOH(B,A,C,N,Fs)

% Simplified by O. David. 

% MVFREQZ multivariate frequency response
% [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF,pCOH2,PDCF,coh,GGC,Af,GPDC] = mvfreqz(B,A,C,N,Fs)
%
% INPUT: 
% ======= 
% A, B	multivariate polynomials defining the transfer function
%
%    a0*Y(n) = b0*X(n) + b1*X(n-1) + ... + bq*X(n-q)
%                          - a1*Y(n-1) - ... - ap*Y(:,n-p)
%
%  A=[a0,a1,a2,...,ap] and B=[b0,b1,b2,...,bq] must be matrices of
%  size  Mx((p+1)*M) and Mx((q+1)*M), respectively. 
%
%  C is the covariance of the input noise X (i.e. D'*D if D is the mixing matrix)
%  N if scalar, N is the number of frequencies 
%    if N is a vector, N are the designated frequencies. 
%  Fs sampling rate [default 2*pi]
% 
%  A,B,C and D can by obtained from a multivariate time series 
%       through the following commands: 
%  [AR,RC,PE] = mvar(Y,P);
%       M = size(AR,1); % number of channels       
%       A = [eye(M),-AR];
%       B = eye(M); 
%       C = PE(:,M*P+1:M*(P+1)); 
%
% OUTPUT: 
% ======= 
% S   	power spectrum
% PDC 	partial directed coherence [2]
% DC  	directed coupling	
% COH 	coherency (complex coherence) [5]
% DTF 	directed transfer function
% pCOH 	partial coherence
% dDTF 	direct Directed Transfer function
% ffDTF full frequency Directed Transfer Function 
% pCOH2 partial coherence - alternative method 
% GGC	a modified version of Geweke's Granger Causality [Geweke 1982]
%	   !!! it uses a Multivariate AR model, and computes the bivariate GGC as in [Bressler et al 2007]. 
%	   This is not the same as using bivariate AR models and GGC as in [Bressler et al 2007]
% Af	Frequency transform of A(z) 
% PDCF 	Partial Directed Coherence Factor [2]
% GPDC 	Generalized Partial Directed Coherence [9,10]
% GdDTF 	Generalized direct Directed Transfer function. Inspired from Baccala et al., 1998

% see also: FREQZ, MVFILTER, MVAR
%
% 
% REFERENCE(S):
% [1] H. Liang et al. Neurocomputing, 32-33, pp.891-896, 2000. 
% [2] L.A. Baccala and K. Samashima, Biol. Cybern. 84,463-474, 2001. 
% [3] A. Korzeniewska, et al. Journal of Neuroscience Methods, 125, 195-207, 2003. 
% [4] Piotr J. Franaszczuk, Ph.D. and Gregory K. Bergey, M.D.
% 	Fast Algorithm for Computation of Partial Coherences From Vector Autoregressive Model Coefficients
%	World Congress 2000, Chicago. 
% [5] Nolte G, Bai O, Wheaton L, Mari Z, Vorbach S, Hallett M.
%	Identifying true brain interaction from EEG data using the imaginary part of coherency.
%	Clin Neurophysiol. 2004 Oct;115(10):2292-307. 
% [6] Schlogl A., Supp G.
%       Analyzing event-related EEG data with multivariate autoregressive parameters.
%       (Eds.) C. Neuper and W. Klimesch, 
%       Progress in Brain Research: Event-related Dynamics of Brain Oscillations. 
%       Analysis of dynamics of brain oscillations: methodological advances. Elsevier. 
% [7] Bressler S.L., Richter C.G., Chen Y., Ding M. (2007)
%	Cortical fuctional network organization from autoregressive modelling of loal field potential oscillations.
%	Statistics in Medicine, doi: 10.1002/sim.2935 
% [8] Geweke J., 1982	
%	J.Am.Stat.Assoc., 77, 304-313.
% [9] L.A. Baccala, D.Y. Takahashi, K. Sameshima. (2006) 
% 	Generalized Partial Directed Coherence. 
%	Submitted to XVI Congresso Brasileiro de Automatica, Salvador, Bahia.  
% [10] L.A. Baccala, D.Y. Takahashi, K. Sameshima. 
% 	Computer Intensive Testing for the Influence Between Time Series, 
%	Eds. B. Schelter, M. Winterhalder, J. Timmer: 
%	Handbook of Time Series Analysis - Recent Theoretical Developments and Applications
%	Wiley, p.413, 2006.

%	$Id: mvfreqz.m 4301 2007-11-26 15:33:54Z schloegl $
%	Copyright (C) 1996-2007 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the TSA-toolbox. See also 
%       http://hci.tugraz.at/schloegl/matlab/tsa/
%       http://octave.sourceforge.net/
%       http://biosig.sourceforge.net/


% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% Boston, MA  02110-1301, USA.

[K1,K2] = size(A);
p = K2/K1-1;
%a=ones(1,p+1);
[K1,K2] = size(B);
q = K2/K1-1;
%b=ones(1,q+1);
if nargin<3
        C = eye(K1,K1);
end;
if nargin<4,
        N = 512;
end;
if nargin<5,
        Fs= 1;        
end;
if all(size(N)==1),	
        f = (0:N-1)*(Fs/(2*N));
else
        f = N;
        N = length(N);
end;
z = i*2*pi/Fs;

h=zeros(K1,K1,N);
S=zeros(K1,K1,N);
COH=zeros(K1,K1,N);

for n=1:N,
        atmp = zeros(K1);
        for k = 1:p+1,
                atmp = atmp + A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;   
                
        btmp = zeros(K1);
        for k = 1:q+1,
                btmp = btmp + B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
        end;        
        h(:,:,n)  = atmp\btmp;        
        S(:,:,n)  = h(:,:,n)*C*h(:,:,n)'/Fs;        
end;

for k1=1:K1;
        for k2=1:K2;
                COH(k1,k2,:)   = (S(k1,k2,:))./sqrt(abs(S(k1,k1,:).*S(k2,k2,:)));
        end;
end;





