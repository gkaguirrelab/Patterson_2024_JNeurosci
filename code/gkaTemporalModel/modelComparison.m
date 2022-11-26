% mergeMRIModelFitResults

% Path to the data dir
dataPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the raw mri data
mriData = loadMRIResponseData();

% Some parameters
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
modelTypes = {'stimulus','cell','mix'};

nBoots = 11;

data = nan(length(subjects),length(modelTypes),nBoots);

for ss = 1:length(subjects)

    % Loop over the model types
    for cc = 1:length(modelTypes)

        % Get the list of result files for this subject
        resultList = dir(fullfile(dataPath,'temporalModelResults',modelTypes{cc},subjects{ss},'mriTemporalModel_*.mat'));

        % Load each result, and save the fVal
        for rr=1:length(resultList)
            loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
            load(loadPath,'mriTemporalModel');
            data(ss,cc,rr)=mriTemporalModel.(subjects{ss}).fVal;
        end

    end

    % Report the comparisons
    [~,p,~,stats] = ttest(squeeze(data(ss,1,:)),squeeze(data(ss,2,:)));
    fprintf([subjects{ss} ': stimulus - cell, t(df) = %2.2f (%d), p = %2.2e \n'],stats.tstat,stats.df,p)

    [~,p,~,stats] = ttest(squeeze(data(ss,1,:)),squeeze(data(ss,3,:)));
    fprintf([subjects{ss} ': stimulus - mix, t(df) = %2.2f (%d), p = %2.2e \n'],stats.tstat,stats.df,p)

end

