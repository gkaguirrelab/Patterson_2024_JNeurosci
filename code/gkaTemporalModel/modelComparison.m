% mergeMRIModelFitResults

% Path to the data dir
dataPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the raw mri data
mriData = loadMRIResponseData();

% Some parameters
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
modelTypes = {'stimulus','cell'};

for ss = 1:length(subjects)

    % Loop over the model types

        % Get the list of result files for this subject
        resultList = dir(fullfile(dataPath,'temporalModelResults','cell',subjects{ss},'mriTemporalModel_*.mat'));

        % Load each result, and save the fVal
        for rr=1:length(resultList)
            % Load the cell result
            loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
            load(loadPath,'mriTemporalModel');
            data(ss,2,rr)=mriTemporalModel.(subjects{ss}).fVal;

            % Load the stimulus result
            loadPath = strrep(loadPath,'/cell/','/stimulus/');
            load(loadPath,'mriTemporalModel');
            data(ss,1,rr)=mriTemporalModel.(subjects{ss}).fVal;
        end


    % Report the comparisons
    [~,p,~,stats] = ttest(squeeze(data(ss,1,:)),squeeze(data(ss,2,:)));
    fprintf([subjects{ss} ': stimulus - cell, t(df) = %2.2f (%d), p = %2.2e \n'],stats.tstat,stats.df,p)


end

