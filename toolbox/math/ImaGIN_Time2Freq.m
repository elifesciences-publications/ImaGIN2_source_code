function freq = ImaGIN_Time2Freq(Time)
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
% Copyright (c) 2000-2017 Inserm
% =============================================================================-
%
% Authors: Olivier David

TimeStep=Time(2)-Time(1);
freqmax	= 1/(2*TimeStep);
if mod(length(Time),2) == 1
    freqpas	= 1/(max(Time)-min(Time));
    freq	= [0:freqpas:freqmax+0.0000001,-freqmax-0.0000001:freqpas:-freqpas];
else
    freqpas	= 1/(max(Time)-min(Time)+TimeStep);
    freq	= [0:freqpas:freqmax+0.0000001,-freqmax-0.0000001+freqpas:freqpas:-freqpas];
    if length(freq) ~= length(Time)
        freq	= [0:freqpas:freqmax+freqpas+0.0000001,-freqmax-0.0000001:freqpas:-freqpas];
    end
end
