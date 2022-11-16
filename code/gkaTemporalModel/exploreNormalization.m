% Where we will save the temporal model results
dataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults');

%% Load the RGC temporal sensitivity model
load(fullfile(dataDir,'rgcTemporalModel.mat'),'rgcTemporalModel');

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

eccDeg = studiedEccentricites(1);
cellClass = 'midget';
stimulusDirection = 'LminusM';

stimulusContrastScale = returnStimulusContrastScale(cellClass,stimulusDirection);
rgcRF=returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,0,stimulusContrastScale,[],[],[],[],[]);

syms f w x
rgIRFSym = ifourier(subs(rgcRF,f,w/(2*pi)),w,x);
deltaT = 0.001;
myTime = 0:deltaT:0.5;
irf = double(subs(rgIRFSym,x,myTime))./1e5;

t1=0.05;
w=0.4;
t2=0.1;
n=1.4;
sigma = 1;

figure
for ff=1:length(studiedFreqs)
    response = stepConvSinStim(irf,deltaT,studiedFreqs(ff));
    response = abs(response);
    responseGamma = stepGammaConv(response,deltaT,t1,0);
%    responseExp = stepExpDecayConv(abs(responseGamma),deltaT,t2);
    responseDenom = sigma^n + abs(responseGamma).^n;
    responseDN = abs(response).^n ./ responseDenom;
    subplot(3,3,ff)
    plot(response,'-k'); hold on
    plot(responseGamma,'-r');plot(responseDN,'-b');
    meanResp(ff)=mean(responseDN);
end

figure
data = mean(mriData.gka.LminusM.ecc1);
semilogx(studiedFreqs,data,'ok');
hold on
semilogx(studiedFreqs,meanResp,'-r');

foo=1;

