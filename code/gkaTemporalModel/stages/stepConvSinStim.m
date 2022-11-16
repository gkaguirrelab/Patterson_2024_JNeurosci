function response = stepConvSinStim(irf,deltaT,freq)

stimDur = 12;
halfCosineRamp = 1.5;

% Make a stimulus
x = 0:deltaT:stimDur-deltaT;
hc = (-cos(x*2*pi*(1/3))+1)/2;
hc(and(x>halfCosineRamp,x<(stimDur-halfCosineRamp)))=1;
stimulus = sin(x*2*pi*freq).*hc;

%padTime = repmat(0,1,round(1/sampleRate));

%stimulus = [stimulus padTime];

% Convolve the stimulus by the irf
temp = conv(stimulus,irf);
response = temp(1:end-length(irf)+1);

end