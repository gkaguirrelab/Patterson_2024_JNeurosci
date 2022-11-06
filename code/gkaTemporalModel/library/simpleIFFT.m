% function [frq, amp, phase] = simpleFFT( signal, ScanRate)
% Purpose: perform an FFT of a real-valued input signal, and generate the single-sided 

% output, in amplitude and phase, scaled to the same units as the input.

%inputs: 

%    signal: the signal to transform

%    ScanRate: the sampling frequency (in Hertz)

% outputs:

%    frq: a vector of frequency points (in Hertz)

%    amp: a vector of amplitudes (same units as the input signal)

%    phase: a vector of phases (in radians)

function [signal, sampleRate] = simpleIFFT( frq, amp, phase)

deltaf = frq(2)-frq(1);

halfn = length(frq);
n = 2*(halfn + 1);

sampleRate = 0.5/frq(end);

theta = [phase phase(end) -fliplr(phase(2:end-1))];
r = [amp(1) amp(2:end-1) amp(end) fliplr(amp(2:end))];

z = r.*exp(1i.*theta);

signal = real(ifft(z)); %do the actual work

end