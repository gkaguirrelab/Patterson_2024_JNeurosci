% Housekeeping
clear

% Model and search type
modelType = 'stimulus';
paramSearch = 'full';

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
subjects = mriFullResultSet.meta.subjects;
plotColor = mriFullResultSet.meta.plotColor;
paramCounts = mriFullResultSet.meta.paramCounts;
cellClasses = {'midget','bistratified','parasol'};
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
nCells = length(cellClasses);
subjectLineSpec = {'-','--'};

% Loop over subjects and calculate the chromaticRG / luminance gain ratio
for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);

    gainVals = {};
    for whichStim = 1:length(stimulusDirections)
        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs + 1;
        gainVals(whichStim) = {mriFullResultSet.(subjects{whichSub}).pMRI(:,startIdx:startIdx+nEccs-1)};
    end

    chromRGToLumRatio(whichSub,:) = mean(gainVals{1}./gainVals{3});
end

subjectGainRatio = mean(chromRGToLumRatio(2,:)./chromRGToLumRatio(1,:));

% Now we will search for an LM ratio for subject 1 that results in a set of
% gain parameters that match subject 2
whichSub = 1;
    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);

    % Set the gain values for chromatic RG to reflect the scaling
    pMRI(5) = pMRI(5)*subjectGainRatio;
    pMRI(17:22) = pMRI(17:22)*subjectGainRatio;

    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1W = 1./std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnW = 1./std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

            results = fitMRIResponse(...
                pMRI,...
                stimulusDirections,studiedEccentricites,studiedFreqs,...
                v1Y,v1W,lgnY,lgnW,...
                'stimulus',false,'LMRatio',true);

