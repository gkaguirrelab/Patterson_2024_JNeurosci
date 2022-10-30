% mergeMRIModelFitResults

% Path to the temporal model results
resultsPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults');

% Loop over subjects
subjects = {'gka','asb'};

% Clear the results structure
mriFullResultSet = [];

for ss = 1:length(subjects)

    % Initialize the fields for this subject
    mriFullResultSet.(subjects{ss}).pMRI = [];
    mriFullResultSet.(subjects{ss}).v1Y = [];
    mriFullResultSet.(subjects{ss}).lgnY = [];

    % Get the list of result files for this subject
    resultList = dir(fullfile(resultsPath,subjects{ss},'mriTemporalModel_*.mat'));

    % Load each result, and merge it into the full set
    for rr=1:length(resultList)
        loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
        load(loadPath,'mriTemporalModel');
        mriFullResultSet.(subjects{ss}).pMRI(end+1,:) = mriTemporalModel.(subjects{ss}).pMRI;
        mriFullResultSet.(subjects{ss}).v1Y(end+1,:) = mriTemporalModel.(subjects{ss}).v1Y;
        mriFullResultSet.(subjects{ss}).lgnY(end+1,:) = mriTemporalModel.(subjects{ss}).lgnY;
    end
end

mriFullResultSet.meta = mriTemporalModel.meta;
savePath = fullfile(resultsPath,'mriFullResultSet.mat');
save(savePath,'mriFullResultSet');