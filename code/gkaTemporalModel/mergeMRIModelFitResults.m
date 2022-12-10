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
paramSearch = 'zeroSurroundIndex';
paramSearch = 'gainOnly';
paramSearch = 'zeroSurroundIndexFreeFilt';
%}
paramSearch = 'full';

% The number of acquisitions obtained for each measurement
nAcqs = 12;
nEcc = 6;

% The number of boot strap resamples of the data to get 95% CIs
nBoots = 10000;

% Index name
indexName = ['mriTemporalModel_' regexprep(num2str(1:12),' +', '-') '_' paramSearch '.mat'];

% Loop over the model types
for cc = 1:length(modelTypes)

    % Clear the results structure
    mriFullResultSet = [];

    for ss = 1:length(subjects)

        % Calculate the CIs on the boot strapped data
        lgnY = []; lgnW = []; v1Y = []; v1W = []; v1AvgY = [];
        for bb=1:10000

            % Get a resample with replacement of the acquisitions
            bootIdx = sort(datasample(1:nAcqs,nAcqs));

            lgnYtemp = []; v1Ytemp = []; v1AvgYtemp = [];
            for whichStim = 1:length(stimulusDirections)

                % Extract the relevant LGN data
                thisMatrix = mriData.(subjects{ss}).(stimulusDirections{whichStim}).lgn(bootIdx,:);
                lgnYtemp = [lgnYtemp mean(thisMatrix)];

                % Extract the relevant V1 data acros eccentricities
                for ee = 1:nEcc
                    thisMatrix = mriData.(subjects{ss}).(stimulusDirections{whichStim}).(['ecc' num2str(ee)])(bootIdx,:);
                    v1Ytemp = [v1Ytemp mean(thisMatrix)];
                    v1EccAccum(ee,:) = mean(thisMatrix);
                end

                % Get the average across eccentricities
                v1AvgYtemp = [v1AvgYtemp, mean(v1EccAccum)];

            end
            lgnY(bb,:) = lgnYtemp;
            v1Y(bb,:) = v1Ytemp;
            v1AvgY(bb,:) = v1AvgYtemp;
        end
        lgnY = sort(lgnY); v1Y = sort(v1Y); v1AvgY = sort(v1AvgY);

        % Initialize the fields for this subject
        mriFullResultSet.(subjects{ss}).pMRI = [];
        mriFullResultSet.(subjects{ss}).v1YMean = mean(v1Y);
        mriFullResultSet.(subjects{ss}).v1W = 1./std(v1Y);
        mriFullResultSet.(subjects{ss}).lgnYMean = mean(lgnY);
        mriFullResultSet.(subjects{ss}).lgnW = 1./std(lgnY);
        mriFullResultSet.(subjects{ss}).v1AvgYMean = mean(v1AvgY);
        mriFullResultSet.(subjects{ss}).v1AvgW = 1./std(v1AvgY);
        mriFullResultSet.(subjects{ss}).v1Y_lowCI = v1Y(round((0.5-0.68/2)*nBoots),:);
        mriFullResultSet.(subjects{ss}).v1AvgY_lowCI = v1AvgY(round((0.5-0.68/2)*nBoots),:);
        mriFullResultSet.(subjects{ss}).lgnY_lowCI = lgnY(round((0.5-0.68/2)*nBoots),:);
        mriFullResultSet.(subjects{ss}).v1Y_highCI = v1Y(round((0.5+0.68/2)*nBoots),:);
        mriFullResultSet.(subjects{ss}).v1AvgY_highCI = v1AvgY(round((0.5+0.68/2)*nBoots),:);
        mriFullResultSet.(subjects{ss}).lgnY_highCI = lgnY(round((0.5+0.68/2)*nBoots),:);

        % Get the list of result files for this subject
        resultList = dir(fullfile(dataPath,'temporalModelResults',modelTypes{cc},subjects{ss},['mriTemporalModel_*_' paramSearch  '.mat']));

        % Load each result, and merge it into the full set
        for rr=1:length(resultList)
            loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
            load(loadPath,'mriTemporalModel');
            mriFullResultSet.(subjects{ss}).pMRI(end+1,:) = mriTemporalModel.(subjects{ss}).pMRI;

            if rr==1 && ss==1
                mriFullResultSet.meta = mriTemporalModel.meta;
            end

            if isfield(mriTemporalModel.(subjects{ss}),'paramCounts')
                mriFullResultSet.meta.modelType = modelTypes{cc};
                mriFullResultSet.meta.paramCounts = mriTemporalModel.(subjects{ss}).paramCounts;
                mriFullResultSet.meta.cellClasses = mriTemporalModel.(subjects{ss}).cellClasses;
            end
        end

    end

    savePath = fullfile(dataPath,'temporalModelResults',modelTypes{cc},['mriFullResultSet_' paramSearch '.mat']);
    save(savePath,'mriFullResultSet');

end
