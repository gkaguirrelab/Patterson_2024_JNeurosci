
% Housekeeping
clear

modelType = 'stimulus';
paramSearch = 'full';

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','v1',modelType);
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

        for whichEcc =1:nEccs+1

            rowIdx = (whichSub-1)*length(stimulusDirections)*(nEccs+1) + (whichStim-1)*(nEccs+1) + whichEcc;

            if whichEcc==(nEccs+1)
                eccName = 'lgn';
            else
                eccName = sprintf('%ddeg',round(studiedEccentricites(whichEcc)));
            end

            T.Properties.RowNames(rowIdx) = {[subjects{whichSub} '_' stimulusDirections{whichStim} '_' eccName]};

            % Unique params
            if whichStim==1 && whichEcc==1
                paramIdx = 1;
                T(rowIdx,1) = {sprintf('%2.2f ± %2.2f',pMRI(paramIdx),pMRISEM(paramIdx))};
                paramIdx = 2;
                T(rowIdx,2) = {sprintf('%2.1f ± %2.1f',pMRI(paramIdx),pMRISEM(paramIdx))};
            else
                T(rowIdx,1:2) = {'',''};
            end

            % Varies by stimulus, fixed across ecentricity
            if whichEcc==1
                paramIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + 1;
                T(rowIdx,3) = {sprintf('%2.1f ± %2.1f',pMRI(paramIdx),pMRISEM(paramIdx))};
            else
                T(rowIdx,3) = {''};
            end

            % Get the weight and gain idx
            if whichEcc==(nEccs+1)
                paramIdxWeight = paramCounts.unique + (whichStim-1)*paramCounts.lgn + 1;
                paramIdxGain = paramCounts.unique + (whichStim-1)*paramCounts.lgn + 2;
            else
                paramIdxWeight = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*0 + whichEcc;
                paramIdxGain = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*1 + whichEcc;

            end
            T(rowIdx,4) = {sprintf('%2.0f ± %2.1f',pMRI(paramIdxWeight)*100,pMRISEM(paramIdxWeight)*100)};
            T(rowIdx,5) = {sprintf('%2.2f ± %2.2f',pMRI(paramIdxGain),pMRISEM(paramIdxGain))};

        end

    end
end

% Save the table
writetable(T,fullfile(savePath,'fullStimulusParamsTable.csv'),'WriteRowNames',true);
