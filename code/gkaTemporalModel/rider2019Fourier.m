% rider2019Fourier
%
% Script to replicate a model and figures from:
%
%   Rider AT, Henning GB, Stockman A. Light adaptation controls visual
%   sensitivity by adjusting the speed and gain of the response to light.
%   PloS one. 2019 Aug 7;14(8):e0220358.
%
% Specifically, the Fourier form of the core model in the appendix (page 2,
% equation "B"), then transformed to the time domain. Compare the output to
% the plots at the bottom right of Figure 8 of the main manuscript.

fc = 15; % corner frequency in Hz for 6 stages
fcl = 30; % corner frequency in Hz for 2 stages
k = 0.8; % relative strength of the "lead compensators" (feedback stages)
g = 10^7.5; % Overall gain

syms f w

% Equation "B" from the appendix of Rider 2019 PLoS One.
mySymFourierF = g * ...         % Overall gain
    (1i.*f+(1-k)*fc).^2 ./ ...  % The "lead compensators"
    ( (1i.*f+fc).^6 .* ...      % 6 low-pass stages with corner fc
    (1i .* f + fcl).^2 );       % 2 low-pass stages with corner fcl

% Convert units from frequency to w
mySymFourierW = subs(mySymFourierF,w / (2*pi));

% Use the inverse fourier transform to obtain the response in time
mySymTime = ifourier(mySymFourierW);

% Show some plots
figure

% The TTF in frequencies
myFreqs = 0:50;
ttfComplex = eval(subs(mySymFourierF,myFreqs));
subplot(3,1,1)
semilogy(myFreqs,abs(ttfComplex),'-r');
axHandle = gca;
ylim([1e-3 1e0]);
axHandle.YTick = [1e-3 1e-2 1e-1 1e0];
xlabel('frequency [Hz]'); ylabel('gain');

% The phase change by frequency, completing the Bode plot
subplot(3,1,2)
plot(myFreqs,angle(ttfComplex)*(180/pi),'*k');
ylim([-180 180]);
axHandle = gca;
axHandle.YTick = [-180 -90 0 90 180];
xlabel('frequency [Hz]'); ylabel('phase [deg]');

% The impulse response function in time
subplot(3,1,3)
myTime = 0:0.001:0.2;
irf = eval(subs(mySymTime,myTime));
plot(myTime*1000,irf,'r')
xlabel('time [msec]'); ylabel('response');

