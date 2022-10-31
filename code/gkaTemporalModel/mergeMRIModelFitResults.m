% mergeMRIModelFitResults

% Path to the data dir
dataPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the raw mri data
mriData = loadMRIResponseData();

% Loop over subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};

% The number of acquisitions obtained for each measurement
nAcqs = 12;
nEcc = 6;

% The number of boot strap resamples of the data to get 95% CIs
nBoots = 10000;

% Clear the results structure
mriFullResultSet = [];

for ss = 1:length(subjects)

    % Calculate the CIs on the boot strapped data
    lgnY = []; lgnW = []; v1Y = []; v1W = [];
    for bb=1:10000

        % Get a resample with replacement of the acquisitions
        bootIdx = sort(datasample(1:nAcqs,nAcqs));

        lgnYtemp = []; v1Ytemp = [];
        for whichStim = 1:length(stimulusDirections)

            % Extract the relevant LGN data
            thisMatrix = mriData.(subjects{ss}).(stimulusDirections{whichStim}).lgn(bootIdx,:);
            lgnYtemp = [lgnYtemp mean(thisMatrix)];

            % Extract the relevant V1 data acros eccentricities
            for ee = 1:nEcc
                thisMatrix = mriData.(subjects{ss}).(stimulusDirections{whichStim}).(['ecc' num2str(ee)])(bootIdx,:);
                v1Ytemp = [v1Ytemp mean(thisMatrix)];
            end
        end
        lgnY(bb,:) = lgnYtemp;
        v1Y(bb,:) = v1Ytemp;
    end

    lgnY = sort(lgnY); v1Y = sort(v1Y);

    % Initialize the fields for this subject
    mriFullResultSet.(subjects{ss}).pMRI = [];
    mriFullResultSet.(subjects{ss}).v1YMean = mean(v1Y);
    mriFullResultSet.(subjects{ss}).lgnYMean = mean(lgnY);
    mriFullResultSet.(subjects{ss}).v1Y_lowCI = v1Y(round((0.5-0.68/2)*nBoots),:);
    mriFullResultSet.(subjects{ss}).lgnY_lowCI = lgnY(round((0.5-0.68/2)*nBoots),:);
    mriFullResultSet.(subjects{ss}).v1Y_highCI = v1Y(round((0.5+0.68/2)*nBoots),:);
    mriFullResultSet.(subjects{ss}).lgnY_highCI = lgnY(round((0.5+0.68/2)*nBoots),:);

    % Get the list of result files for this subject
    resultList = dir(fullfile(dataPath,'temporalModelResults',subjects{ss},'mriTemporalModel_*.mat'));

    % Load each result, and merge it into the full set
    for rr=1:length(resultList)
        loadPath = fullfile(resultList(rr).folder,resultList(rr).name);
        load(loadPath,'mriTemporalModel');
        mriFullResultSet.(subjects{ss}).pMRI(end+1,:) = mriTemporalModel.(subjects{ss}).pMRI;
    end
end

mriFullResultSet.meta = mriTemporalModel.meta;
savePath = fullfile(dataPath,'temporalModelResults','mriFullResultSet.mat');
save(savePath,'mriFullResultSet');