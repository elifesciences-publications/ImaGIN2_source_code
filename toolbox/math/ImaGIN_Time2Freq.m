function freq = ImaGIN_Time2Freq(Time)

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
