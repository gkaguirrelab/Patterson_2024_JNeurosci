% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','mriTemporalModel.mat');
load(loadPath,'mriTemporalModel');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriTemporalModel.meta.studiedFreqs;
studiedEccentricites = mriTemporalModel.meta.studiedEccentricites;
subjects = mriTemporalModel.meta.subjects;
stimulusDirections = mriTemporalModel.meta.stimulusDirections;
plotColor = mriTemporalModel.meta.plotColor;
nFixedParams = mriTemporalModel.meta.nFixedParams;
nFloatByEccParams = mriTemporalModel.meta.nFloatByEccParams;
nUniqueParams = mriTemporalModel.meta.nUniqueParams;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
cellClassOrder = {'midgetChrom','parasol','bistratified','midgetAchrom'};

for whichSub = 1:length(subjects)

    pMRI = mriTemporalModel.(subjects{whichSub}).pMRI;
    v1Y = mriTemporalModel.(subjects{whichSub}).data.v1Y;
    lgnY = mriTemporalModel.(subjects{whichSub}).data.lgnY;

    [~,v1YFitMatrix] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,nUniqueParams,nFixedParams);
    lgnYFit = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,nUniqueParams,nFixedParams);

    figure

    for whichStim = 1:length(stimulusDirections)


%        Show the cortical filter
        secondOrderFc = pMRI( (whichStim-1)*(nUniqueParams/3)+1 );
        secondOrderQ = pMRI( (whichStim-1)*(nUniqueParams/3)+2 );
        syms f; rf = stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
        myFreqs = logspace(log10(0.5),log10(100),101);
        ttfComplex = double(subs(rf,myFreqs));
        gainVals = abs(ttfComplex);
        semilogx(myFreqs,gainVals,['-' plotColor{whichStim}]);
        hold on
        xlabel('frequency [Hz]'); ylabel('gain');


    end

end
