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
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults');
load(fullfile(loadPath,'mriFullResultSet.mat'),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
nFixedParams = mriFullResultSet.meta.nFixedParams;
nFloatByEccParams = mriFullResultSet.meta.nFloatByEccParams;
nUniqueParams = mriFullResultSet.meta.nUniqueParams;
postReceptoralPaths = mriFullResultSet.meta.postReceptoralPaths;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
chromAchromIndex = [1 1 2];
nParamsPerCellBlock = nFixedParams+nEccs*2;

eccDeg = 20.5;
ee = 4;

whichSub = 1;

pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);

% Assemble the LGN reponse params
lgnSurroundDelay = pMRI(1);
lgnSurroundIndex = pMRI(2);

% Loop over the stims
for ss=1:length(stimulusDirections)

    figHandleLGN{ss} = figure('Position',  [100, 100, 200, 400]);
    figHandleV1{ss} = figure('Position',  [100, 100, 200, 400]);

    % Identify which cell classes are relevant for this stimulus direction
    switch stimulusDirections{ss}
        case 'LminusM'
            cellClasses = {'midget'};
            plotColor = {'-r'};
            pathwayIndex = 1;
            chromAchromIndex = 1;
        case 'S'
            cellClasses = {'bistratified'};
            plotColor = {'-b'};
            pathwayIndex = 2;
            chromAchromIndex = 1;
        case 'LMS'
            cellClasses = {'parasol','midget'};
            plotColor = {'-k','-.k'};
            pathwayIndex = [3 4];
            chromAchromIndex = 2;
    end

    for cc = 1:length(cellClasses)

        % Grab the LGN parameters, which are organized by RGC class
        switch cellClasses{cc}
            case 'midget'
                lgnGain = pMRI(3);
            case 'bistratified'
                lgnGain = pMRI(4);
            case 'parasol'
                lgnGain = pMRI(5);
        end

        % Grab the V1 MRI parameters that are fixed across
        % eccentricity, but vary with stimulus class (chromatic or
        % achromatic)
        secondOrderFc = pMRI(5+(chromAchromIndex-1)*2+1);
        secondOrderQ = pMRI(5+(chromAchromIndex-1)*2+2);

        % Grab the V1 MRI parameters that are fixed across
        % eccentricity, but vary with post-receptoral path
        v1SurroundDelay = pMRI(nUniqueParams+(pathwayIndex(cc)-1)*nParamsPerCellBlock+1);

        % Grab the "floating" parameters that vary with eccentricity
        v1SurroundIndex = pMRI(nUniqueParams+(pathwayIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+ee);
        v1Gain = pMRI(nUniqueParams+(pathwayIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+nEccs+ee);

        % Assemble the staged parameters
        surroundDelay = [lgnSurroundDelay v1SurroundDelay];
        surroundIndex = [lgnSurroundIndex v1SurroundIndex];
        gain = [lgnGain v1Gain];

        rfPostRetinalLGN{ss,cc} = ...
            returnPostRetinalRF(cellClasses{cc},stimulusDirections{ss},...
            rgcTemporalModel,eccDeg,1,...
            lgnSurroundDelay,lgnSurroundIndex,lgnGain,[],[]);

        rfPostRetinalV1{ss,cc} = ...
            returnPostRetinalRF(cellClasses{cc},stimulusDirections{ss},...
            rgcTemporalModel,eccDeg,2,...
            surroundDelay,surroundIndex,gain,secondOrderFc,secondOrderQ);

        plotRF(rfPostRetinalLGN{ss,cc},figHandleLGN{ss},plotColor{cc});
        plotRF(rfPostRetinalV1{ss,cc},figHandleV1{ss},plotColor{cc});
    end

    plotName = ['lgnRF_' subjects{1} '_' stimulusDirections{ss} '_' num2str(eccDeg,2) '_ModelFit.pdf' ];
    saveas(figHandleLGN{ss},fullfile(savePath,plotName));

    plotName = ['V1RF_' subjects{1} '_' stimulusDirections{ss} '_' num2str(eccDeg,2) '_ModelFit.pdf' ];
    saveas(figHandleV1{ss},fullfile(savePath,plotName));

end