% Where we will save the temporal model results
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults');


%% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);


%% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEcc = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k','g'};
postReceptoralPaths = {'midget.LminusM','bistratified.S','parasol.LMS','midget.LMS'};

% The number of acquisitions obtained for each measurement
nAcqs = 12;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

eccDeg = studiedEccentricites(6);
cellClass = 'midget';
stimulusDirection = 'LminusM';

[rfPostRetinal, rfRGC] = returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,1,0,0,5e-2,[],[]);

myFreqs = linspace(0,1000,401);
ttfComplex = double(subs(rfPostRetinal,myFreqs));
[irf, sampleRate] = simpleIFFT( myFreqs, abs(ttfComplex), angle(ttfComplex));
myTime = 0:sampleRate:(length(irf)-1)*sampleRate;
irf = irf-irf(1);

% Apply a gamma kernel to the irf
gammaTime = 0:sampleRate:1;
t1 = 0.025; % Midget LMS ecc 6  = 0.025
%w = 0.4;
gammaKernel = gammaTime.*exp(-gammaTime/t1);% - w*gammaTime.*exp(-gammaTime/(1.5*t1)) ;
gammaKernel = gammaKernel / sum(gammaKernel);
temp = conv(irf,gammaKernel);
irfGamma = temp(1:end-length(gammaKernel)+1);
irfGamma = irfGamma/sum(irfGamma);

% Loop over frequencies
fitFrequencies = logspace(log10(2),log10(64),20);
response = [];

for ff = 1:length(fitFrequencies)

%% Make a stimulus

x = 0:sampleRate:12-sampleRate;
padTime = repmat(0,1,round(1/sampleRate));
hc = (-cos(x*2*pi*(1/3))+1)/2;
hc(and(x>1.5,x<10.5))=1;
stimulus = sin(x*2*pi*fitFrequencies(ff)).*hc;
stimulus = [stimulus padTime];

%% Convolve the stimulus by the irf
rgcResponse = conv(stimulus,irfGamma);
rgcResponse = rgcResponse(1:end-length(irfGamma)+1);

n=2; % Midget LMS ecc 6  = 1.4
t2 = 0.03; % Midget LMS ecc 6  = 0.07 
sigma = 0.01; % Midget LMS ecc 6  =  0.4
expTime = 0:sampleRate:5;
expKernel = exp(-expTime/t2);
expStart = find(irf>max(irf)/1e3,1);
expKernel = circshift(expKernel,expStart);
expKernel(1:expStart)=0;
expKernel = expKernel./sum(expKernel);
temp = conv(rgcResponse,expKernel);
denom = temp(1:end-length(expKernel)+1);

temporalResponse = abs(rgcResponse).^n ./ (sigma.^n + abs(denom).^n);
response(ff) = mean(temporalResponse);

end

data = mean(mriData.gka.LminusM.ecc6);
semilogx(studiedFreqs,data,'or');
hold on
scaledResponse = (response/max(response))*max(data);
semilogx(fitFrequencies,scaledResponse,'-r');