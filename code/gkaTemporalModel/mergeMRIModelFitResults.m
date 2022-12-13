% mergeMRIModelFitResults

% Path to the data dir
dataPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the raw mri data
mriData = loadMRIResponseData();

% Some parameters
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
modelTypes = {'stimulus'};
%{
paramSearch = 'full';
paramSearch = 'gainOnly';
paramSearch = 'noSurround';
%}
paramSearch = 'full';

% The number of acquisitions obtained for each measurement
nAcqs = 12;
nEcc = 6;

% Index name
indexName = ['mriTemporalModel_' regexprep(num2str(1:12),' +', '-') '_' paramSearch '.mat'];

% Loop over the model types
for cc = 1:length(modelTypes)

    % Clear the results structure
    mriFullResultSet = [];

    for ss = 1:length(subjects)

        % Get the list of result files for this subject
        resultList = dir(fullfile(dataPath,'temporalModelResults',modelTypes{cc},subjects{ss},['mriTemporalModel_*_' paramSearch  '.mat']));

        % Load each result, and merge it into the full set
        for rr=1:length(resultList)
            loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
            load(loadPath,'mriTemporalModel');

            % Initialize the container for the assembled results
            if rr==1
                mriFullResultSet.meta = mriTemporalModel.meta;
                mriFullResultSet.meta.modelType = modelTypes{cc};
                mriFullResultSet.meta.paramCounts = mriTemporalModel.(subjects{ss}).paramCounts;
                mriFullResultSet.meta.cellClasses = mriTemporalModel.(subjects{ss}).cellClasses;
                mriFullResultSet.(subjects{ss}).pMRI = [];
                mriFullResultSet.(subjects{ss}).v1Y = [];
                mriFullResultSet.(subjects{ss}).lgnY = [];
            end

            mriFullResultSet.(subjects{ss}).pMRI(end+1,:) = mriTemporalModel.(subjects{ss}).pMRI;
            mriFullResultSet.(subjects{ss}).v1Y(end+1,:) = mriTemporalModel.(subjects{ss}).v1Y;
            mriFullResultSet.(subjects{ss}).lgnY(end+1,:) = mriTemporalModel.(subjects{ss}).lgnY;

        end

    end

    mriFullResultSet.(subjects{ss}).v1YSEM = std(mriFullResultSet.(subjects{ss}).v1Y);
    mriFullResultSet.(subjects{ss}).lgnYSEM = std(mriFullResultSet.(subjects{ss}).lgnY);

    savePath = fullfile(dataPath,'temporalModelResults',modelTypes{cc},['mriFullResultSet_' paramSearch '.mat']);
    save(savePath,'mriFullResultSet');

end
