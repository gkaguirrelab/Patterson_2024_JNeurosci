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

p0 = [0.01,1.4,0.07,0.4];

data = mean(mriData.gka.LminusM.ecc6);
myResponse = @(p) returnPostRetinalIRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,p);
myObj = @(p) norm(data-myResponse(p));

myResponse(p0);

p=fmincon(myObj,p0);

semilogx(studiedFreqs,data,'or');
hold on
scaledResponse = (response/max(response))*max(data);
semilogx(fitFrequencies,scaledResponse,'-r');