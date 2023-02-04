% mergeMRIModelFitResults

% Path to the data dir
dataPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data');
dataPath = fullfile(dataPath,'temporalModelResults','v1','crossVal');

% Some parameters
subjects = {'gka','asb'};
modelTypes = {'stimulus','cell'};
paramTypes = {'gainOnly','full'};
fieldType = {'lgnR2','v1R2','v1EccR2'};

data = [];

% Get the list of result files for this subject
resultList = dir(fullfile(dataPath,'crossValResults_*.mat'));

% Load each result, and save the fVal
for rr=1:length(resultList)

    % Load the cell result
    loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
    load(loadPath,'crossValResults');

    for ss=1:length(subjects)
        for mm=1:length(modelTypes)
            for pp=1:length(paramTypes)
                for ff=1:length(fieldType)
                    data.(fieldType{ff}).(subjects{ss}).(paramTypes{pp}).(modelTypes{mm})(rr) = ...
                        crossValResults.(subjects{ss}).(modelTypes{mm}).(paramTypes{pp}).(fieldType{ff});
                end
            end
        end
    end
end

%
%     % Report the comparisons
%     [~,p,~,stats] = ttest(squeeze(data(ss,1,:)),squeeze(data(ss,2,:)));
%     fprintf([subjects{ss} ': stimulus - cell, t(df) = %2.2f (%d), p = %2.2e \n'],stats.tstat,stats.df,p)
%
%     nStimBetter = sum(squeeze(data(ss,1,:)) < squeeze(data(ss,2,:)));
%     percentStimBetter = 100*mean((squeeze(data(ss,2,:)) - squeeze(data(ss,1,:)) ) ./ squeeze(data(ss,2,:)));
%     fprintf([subjects{ss} ': stimulus better then cell: %d / %d cases\n'],nStimBetter,length(resultList))
%     fprintf([subjects{ss} ': stimulus %2.0f%% better fit than cell\n'],percentStimBetter)



