function response = returnPostRetinalIRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,p)

% Unpack the parameters
t1 = p(1);
n = p(2); % Midget LMS ecc 6  = 1.4
t2 = p(3); % Midget LMS ecc 6  = 0.07 
sigma = p(4); % Midget LMS ecc 6  =  0.4

% Obtain the chromatic weights
[chromaticCenterWeight,chromaticSurroundWeight] = ...
    returnRGCChromaticWeights(cellClass,stimulusDirection,eccDeg,rgcTemporalModel.LMRatio);

% Obtain the pRGC model params for this cell class and eccentricity
cellIndex = strcmp(cellClass,{'midget','parasol','bistratified'});
pBlockRGC = nan(1,7);
for ii = 1:7
    pBlockRGC(ii) = rgcTemporalModel.pFitByEccen{ii,cellIndex}(eccDeg);
end

% Get the temporal receptive field for this RGC
[rfRGC, rfBipolar, rfCone] = ...
    returnRGCRF(pBlockRGC,rgcTemporalModel.cfCone,rgcTemporalModel.coneDelay,chromaticCenterWeight,chromaticSurroundWeight);

% Convert from response units of gain to impulses/sec by accounting for the
% contrast of the stimulus, which is hard-coded in this called function
rfRGC = rfRGC*returnStimulusContrastScale(cellClass,stimulusDirection);

% Copy the RGC model into the post-retinal variables
rfPostRetinal = rfRGC;

% Drasdo 2007 equation for the midget fraction as a function of
% eccentricity
midgetFraction = @(eccDeg) 0.8928*(1+eccDeg/41.03).^(-1);

% Bistratified RGCs as a function of eccentricity. Values taken from Figure
% 13B of:
%   Dacey DM. Morphology of a small-field bistratified ganglion cell type
%   in the macaque and human retina. Visual neuroscience. 1993
%   Nov;10(6):1081-98.
% Further, set the bistratified density to 0 below 0.5 degrees to model the
% tritanopic foveola.
%{
    daceyXmm = [1.0082    1.9945    2.9589    3.9671    4.9534    5.9616    6.9699    7.9123    8.9205    9.9507   10.9589   11.9671    12.9973   13.9616   14.9699];
    daceyY = 0.01 .* [1.3912    1.5573    2.0519    2.7377    3.2223    3.4308    3.5179    3.7927    4.1928    4.7528    5.3876    5.8085  5.8085    6.1072    6.7515];
    daceyXdeg = 0.1 + 3.4.*daceyXmm + 0.035.*daceyXmm.^2;
    pp = polyfit(daceyXdeg,daceyY,2)
    bistratifiedFraction = @(eccDeg) (eccDeg>0.5).*( pp(1).*eccDeg.^2+pp(2).*eccDeg+pp(3));
    figure
    semilogy(daceyXdeg,daceyY,'ok'); hold on
    semilogy(0:1:90,bistratifiedFraction(0:1:90),'-r')
%}
pp(1) = -0.000002818937309; pp(2) = 0.001127464065019; pp(3) = 0.009832178782458;
bistratifiedFraction = @(eccDeg) (eccDeg>0.5).*( pp(1).*eccDeg.^2+pp(2).*eccDeg+pp(3) );

% At each eccentricity, what fraction of the total number of RGCs is of a
% given cell class? The parasol denisty is what is left over after we
% account for midget and bistratified classes.
switch cellClass
    case 'midget'
        proportionFunc = @(ecc) midgetFraction(ecc);
    case 'parasol'
        proportionFunc = @(ecc) 1.0-midgetFraction(ecc)-bistratifiedFraction(ecc);
    case 'bistratified'
        proportionFunc = @(ecc) bistratifiedFraction(ecc);
end

% Loop over the meridians and obtain the RGCf density functions
meridianAngles = [0 90 180 270];
for mm = 1:4
    totalRGCfDensity{mm} = watsonTotalRFDensityByEccenDegVisual(meridianAngles(mm));
end

% Obtain the mean, total ganglion receptive field density as a function of
% eccentricity
totalRGCfDensityAtEcc = @(eccDeg) mean(cellfun(@(x) x(eccDeg),totalRGCfDensity));

% Scale by cell density at this eccentricity, accounting for the
% change in surface area as a function of eccentricity (eccDeg/90), and
% the fraction of cells of each class.
rfPostRetinal = rfPostRetinal ...
    .* totalRGCfDensityAtEcc(eccDeg) .* (eccDeg/90) ...
    .* proportionFunc(eccDeg);

% Obtain the IRF for this temporal receptive field
myFreqs = linspace(0,1000,401);
ttfComplex = double(subs(rfPostRetinal,myFreqs));
[postRetinalIRF, sampleRate] = simpleIFFT( myFreqs, abs(ttfComplex), angle(ttfComplex));
postRetinalIRF = postRetinalIRF-postRetinalIRF(1);

% Apply a gamma kernel to the irf
gammaTime = 0:sampleRate:1;
gammaKernel = gammaTime.*exp(-gammaTime/t1);
gammaKernel = gammaKernel / sum(gammaKernel);
temp = conv(postRetinalIRF,gammaKernel);
postRetinalIRF = temp(1:end-length(gammaKernel)+1);

% Loop over frequencies
fitFrequencies = [2 4 8 16 32 64];%logspace(log10(2),log10(64),20);
response = [];

for ff = 1:length(fitFrequencies)

% Make a stimulus
x = 0:sampleRate:12-sampleRate;
padTime = repmat(0,1,round(1/sampleRate));
hc = (-cos(x*2*pi*(1/3))+1)/2;
hc(and(x>1.5,x<10.5))=1;
stimulus = sin(x*2*pi*fitFrequencies(ff)).*hc;
stimulus = [stimulus padTime];

% Convolve the stimulus by the irf
rgcResponse = conv(stimulus,postRetinalIRF);
rgcResponse = rgcResponse(1:end-length(postRetinalIRF)+1);

expTime = 0:sampleRate:5;
expKernel = exp(-expTime/t2);
expStart = find(postRetinalIRF>max(postRetinalIRF)/1e3,1);
expKernel = circshift(expKernel,expStart);
expKernel(1:expStart)=0;
expKernel = expKernel./sum(expKernel);
temp = conv(rgcResponse,expKernel);
denom = temp(1:end-length(expKernel)+1);

temporalResponse = abs(rgcResponse).^n ./ (sigma.^n + abs(denom).^n);
response(ff) = mean(temporalResponse);

end

end
