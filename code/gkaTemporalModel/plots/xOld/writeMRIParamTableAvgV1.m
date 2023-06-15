
% Housekeeping
clear

modelType = 'stimulus';
paramSearch = 'full';

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_FullModel',modelType);

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
paramCounts = mriFullResultSet.meta.paramCounts;
cellClasses = {'midget','bistratified','parasol'};
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
fitFrequencies = logspace(0,2,50);
nFreqsForPlotting = length(fitFrequencies);
nCells = length(cellClasses);
subjectLineSpec = {'-','--'};


varNames = {'n','fc [Hz]','d [msecs]','w [%]','g [x10^4]'};
T = table('Size',[2*3*7,length(varNames)],'VariableTypes',repmat({'string'},1,length(varNames)));
T.Properties.VariableNames=varNames;

for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);

    for whichStim = 1:length(stimulusDirections)

        rowIdx = (whichSub-1)*length(stimulusDirections)*(2) + (whichStim-1)*(2) + 1;
        eccName = 'lgn';
        T.Properties.RowNames(rowIdx) = {[subjects{whichSub} '_' stimulusDirections{whichStim} '_' eccName]};

        paramIdxWeight = []; paramIdxGain = [];
        paramIdx = 1;
        T(rowIdx,1) = {sprintf('%2.2f ± %2.2f',pMRI(paramIdx),pMRISEM(paramIdx))};
        paramIdx = 2;
        T(rowIdx,2) = {sprintf('%2.1f ± %2.1f',pMRI(paramIdx),pMRISEM(paramIdx))};
        paramIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + 1;
        T(rowIdx,3) = {sprintf('%2.1f ± %2.1f',pMRI(paramIdx),pMRISEM(paramIdx))};
        paramIdxWeight = paramCounts.unique + (whichStim-1)*paramCounts.lgn + 1;
        paramIdxGain = paramCounts.unique + (whichStim-1)*paramCounts.lgn + 2;
        T(rowIdx,4) = {sprintf('%2.0f ± %2.1f',pMRI(paramIdxWeight)*100,pMRISEM(paramIdxWeight)*100)};
        T(rowIdx,5) = {sprintf('%2.2f ± %2.2f',pMRI(paramIdxGain),pMRISEM(paramIdxGain))};

        rowIdx = (whichSub-1)*length(stimulusDirections)*(2) + (whichStim-1)*(2) + 2;
        eccName = 'avgV1';
        T.Properties.RowNames(rowIdx) = {[subjects{whichSub} '_' stimulusDirections{whichStim} '_' eccName]};
        for rr=1:5
            T(rowIdx,rr) = {' - '};
        end

        paramIdxWeight = []; paramIdxGain = [];
        for whichEcc =1:nEccs
            paramIdxWeight(whichEcc) = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*0 + whichEcc;
            paramIdxGain(whichEcc) = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*1 + whichEcc;
        end
        T(rowIdx,4) = {sprintf('%2.0f ± %2.1f',mean(pMRI(paramIdxWeight))*100,mean(pMRISEM(paramIdxWeight))*100)};
        T(rowIdx,5) = {sprintf('%2.2f ± %2.2f',mean(pMRI(paramIdxGain)),mean(pMRISEM(paramIdxGain)))};

    end

end

% Save the table
writetable(T,fullfile(savePath,'fullStimulusParamsTableAvgV1.csv'),'WriteRowNames',true);
