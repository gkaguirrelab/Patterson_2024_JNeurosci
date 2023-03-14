
%% Housekeeping
clear
close all


%% Initialize the random seed so that we get different perms each run
rng('shuffle');

% Which cortical region to fit
corticalRegion = 'v1';

% What model type do we want? By cell or by stimulus?
modelTypes = {'stimulus'};

% Which set of parameters will we investigate in the bootstrap analysis?
paramSearch = 'fullCortex';

% How many bootstrap resamplings of the data to conduct
nBoots = 1;

% Verbose?
verbose = true;

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
plotColor = {'r','b','k'};

% The number of acquisitions obtained for each measurement
nAcqs = 12;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Define the temporal model that will hold the results
mriTemporalModel.meta.studiedFreqs = studiedFreqs;
mriTemporalModel.meta.studiedEccentricites = studiedEccentricites;
mriTemporalModel.meta.subjects = subjects;
mriTemporalModel.meta.stimulusDirections = stimulusDirections;
mriTemporalModel.meta.plotColor = plotColor;
mriTemporalModel.meta.nBoots = nBoots;
mriTemporalModel.meta.paramSearch = paramSearch;

% Loop over subjects
for whichSub = [1 2]

    % Assemble the data
    lgnY = []; lgnW = []; cortexY = []; cortexW = [];
    for whichStim = 1:length(stimulusDirections)

        % Extract the V1 response across eccentricities
        for ee = 1:nEccs
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).([corticalRegion '_ecc' num2str(ee)]);
            cortexY = [cortexY mean(thisMatrix)];
            cortexW = [cortexW 1./std(thisMatrix)];
        end
    end % loop over stims

    % Loop over models
    for mm=1:length(modelTypes)

        % Load a search seed
        pMRI0 = storedSearchSeeds(whichSub,modelTypes{mm},'full');

        % Report our progress
        curTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
        str=[curTime ' - Subject: ' subjects{whichSub} ', ' modelTypes{mm}];
        fprintf(str);

        % Let's see how long this takes
        tic

        % BADS it.
        results = fitMRIResponse(...
            pMRI0,...
            stimulusDirections,studiedEccentricites,studiedFreqs,...
            cortexY,cortexW,lgnY,lgnW,...
            modelTypes{mm},false,paramSearch,verbose);
        results.pMRI0 = pMRI0;

        % Report our search time and outcome
        searchTimeSecs = toc();
        str=[sprintf('fVal = %2.2f, search time (mins) = %2.1f',results.fVal,searchTimeSecs/60)  '\n'];
        fprintf(str);

        % Report the parameters
        paramCounts = results.paramCounts;
        if verbose
            str = 'pMRI0 = [ ...\n';
            str = [str sprintf(repmat('%2.10f, ',1,paramCounts.unique),results.pMRI(1:paramCounts.unique)) ' ...\n'];
            str = [str sprintf(repmat('%2.10f, ',1,paramCounts.lgn*3),results.pMRI(paramCounts.unique+1:paramCounts.unique+paramCounts.lgn*3)) ' ...\n'];
            for ss=1:length(stimulusDirections)
                startIdx = paramCounts.unique+ paramCounts.lgn*3 + (ss-1)*paramCounts.v1total;
                str = [str sprintf(repmat('%2.10f, ',1,paramCounts.v1total),results.pMRI(startIdx+1:startIdx+paramCounts.v1total)) ' ...\n'];
            end
            str = [str ']; \n'];
            fprintf(str);
        end

    end % model types

end % subjects

