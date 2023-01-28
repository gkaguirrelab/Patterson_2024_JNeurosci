
%% Housekeeping
clear
close all


%% Initialize the random seed so that we get different perms each run
rng('shuffle');


% Some model settings
useMonotonicConstraint = false;
corticalRegion = 'v1';
modelType = {'cell','stimulus'};

% The types of parameters to use
paramSearch = {'gainOnly','full'};

% How many bootstrap cross-validations of the data to conduct
nSplits = 25;

% Verbose?
verbose = false;

% Where we will save the cross validation results?
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults',corticalRegion);

% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};

% The number of acquisitions obtained for each measurement
nAcqs = 12;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Define the structure that will hold the results
crossValResults_fVal.meta.nSplits = nSplits;

% Loop over bootstraps
for bb = 1:nSplits

    splitVec = randperm(12);

    % Loop over subjects
    for whichSub = [1 2]

        % Loop over model type
        for mm = 1:length(modelType)

            % Loop over param searches
            for pp = 1:length(paramSearch)

                % Report our progress
                curTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
                str=[curTime ' - Subject: ' subjects{whichSub} ', ' modelType{mm} ', ' paramSearch{pp} sprintf(', split: %d ...',bb)];
                fprintf(str);

                % Let's see how long this takes
                tic

                % Clear the pMRI0
                pMRI0 = [];

                % Loop over train and test
                for tt=1:2

                    % Get the train or test half of the data
                    bootIdx = splitVec((tt-1)*6+1:tt*6);

                    % Assemble the data
                    lgnY = []; lgnW = []; cortexY = []; cortexW = [];
                    for whichStim = 1:length(stimulusDirections)

                        % Extract the  LGN data
                        thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).lgn(bootIdx,:);
                        lgnY = [lgnY mean(thisMatrix)];

                        % Extract the V1 response across eccentricities
                        for ee = 1:nEccs
                            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).([corticalRegion '_ecc' num2str(ee)])(bootIdx,:);
                            cortexY = [cortexY mean(thisMatrix)];
                        end
                    end

                    % Set the weights equal to unity
                    lgnW = ones(size(lgnY));
                    cortexW = ones(size(cortexY));

                    % Load a search seed
                    switch tt
                        case 1
                            pMRI0 = storedSearchSeeds(whichSub,modelType{mm},paramSearch{pp});
                            thisParamSearch = paramSearch{pp};
                        case 2
                            pMRI0 = pMRI0; % The solution on the last loop
                            thisParamSearch = 'lockAll';
                    end

                    % BADS it.
                    results = fitMRIResponse(...
                        pMRI0,...
                        stimulusDirections,studiedEccentricites,studiedFreqs,...
                        cortexY,cortexW,lgnY,lgnW,...
                        modelType{mm},useMonotonicConstraint,thisParamSearch,verbose);
                    pMRI0 = results.pMRI;

                end

                % Report our search time and outcome
                searchTimeSecs = toc();
                str=[sprintf('Cross-val lgnR2 = %2.2f, v1R2 = %2.2f, search time (mins) = %2.1f',results.lgnR2,results.v1R2,searchTimeSecs/60)  '\n'];
                fprintf(str);

                % Store the value
                crossValResults_fVal.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp})(bb) = results.fVal;
                crossValResults_v1R2.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp})(bb) = results.v1R2;
                crossValResults_lgnR2.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp})(bb) = results.lgnR2;

                % Save after each iteration
                saveSpot = fullfile(saveDir,'crossValResults.mat');
                save(saveSpot,'crossValResults_fVal','crossValResults_v1R2','crossValResults_lgnR2');

            end % paramSearch types

        end % model types

    end % subjects

end % splits


