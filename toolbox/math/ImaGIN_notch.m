function x = ImaGIN_notch(x, varargin)
% USAGE:  x = ImaGIN_notch(x, Fs, Fp, Fmax)
%         x = ImaGIN_notch(x, Wo, BW)

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
%          Francois Tadel, 2017 (code reformatting)

% CALL:  x = ImaGIN_notch(x,Fs,Fp,Fmax)
if (nargin == 4)
    Fs   = varargin{1};
    Fp   = varargin{2};
    Fmax = varargin{3};
    
    Fmax = min([Fmax Fs/2]);
    for i1 = 1:floor(Fmax/Fp)
        Wo = Fp*2/Fs; 
        BW = Wo/35;
        [bnotch,anotch] = local_notch(Wo,BW);
        x = filter(bnotch,anotch,x);
    end
    
% CALL:  x = ImaGIN_notch(x,Wo,BW)
elseif (nargin == 3)
    Wo = varargin{1};
    BW = varargin{2};
    
    [bnotch,anotch] = local_notch(Wo,BW);
    x = filter(bnotch,anotch,x);
end
end


%% ===== FILTER FUNCTION =====
function [num, den] = local_notch(Wo, BW)
    BW = BW * pi;
    Wo = Wo * pi;

    Gb   = 10^(-abs(10*log10(.5))/20);
    beta = (sqrt(1-Gb.^2)/Gb) * tan(BW/2);
    gain = 1/(1+beta);

    num  = gain * [1, -2*cos(Wo), 1];
    den  = [1, -2*gain*cos(Wo), (2*gain-1)];
end


