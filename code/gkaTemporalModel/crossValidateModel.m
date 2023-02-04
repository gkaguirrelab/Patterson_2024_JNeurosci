
%% Housekeeping
clear
close all


%% Initialize the random seed so that we get different perms each run
rng('shuffle');


% Some model settings
useMonotonicConstraint = false;
corticalRegion = 'v1';

% Which model types to test
modelType = {'cell','stimulus'};

% The types of parameters to use
paramSearch = {'gainOnly','full'};

% How many bootstrap cross-validations of the data to conduct
nSplits = 6;

% Verbose?
verbose = false;

% Where we will save the cross validation results?
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults',corticalRegion,'crossVal');

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
crossValResults.meta.nSplits = nSplits;

% Open the parpool
gcp;

% Loop over bootstraps
for bb = 1:nSplits

    % Clear the results structure
    crossValResults = [];
    
    % Define the fit and test sets
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
                        lgnW = [lgnW 1./std(thisMatrix)];

                        % Extract the V1 response across eccentricities
                        for ee = 1:nEccs
                            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).([corticalRegion '_ecc' num2str(ee)])(bootIdx,:);
                            cortexY = [cortexY mean(thisMatrix)];
                            cortexW = [cortexW 1./std(thisMatrix)];
                        end
                    end

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
                str=[sprintf('Cross-val lgnR2 = %2.2f, v1R2 = %2.2f, v1EccR2 = %2.2f, search time (mins) = %2.1f',results.lgnR2,results.v1R2,results.v1EccR2,searchTimeSecs/60)  '\n'];
                fprintf(str);

                % Store the value
                crossValResults.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp}).fVal = results.fVal;
                crossValResults.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp}).lgnR2 = results.lgnR2;
                crossValResults.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp}).v1R2 = results.v1R2;
                crossValResults.(subjects{whichSub}).(modelType{mm}).(paramSearch{pp}).v1EccR2 = results.v1EccR2;

            end % paramSearch types

        end % model types

    end % subjects

    % Save after each iteration
    bootLabel = regexprep(num2str(splitVec),' +', '-');
    saveSpot = fullfile(saveDir,['crossValResults_' bootLabel '.mat']);
    save(saveSpot,'crossValResults');

    % After every loop, close the parpool and re-open. This is a kludge to
    % handle a quasi-memory link that occurs in the muPad symbolic toolbox,
    % which retains a copy of every symbolic variable, causing a growing
    % memory overhead. We also clear out old / bad jobs from the cluster.
    delete(gcp('nocreate'));
    myCluster = parcluster('Processes');
    delete(myCluster.Jobs)
    gcp;

end % splits


