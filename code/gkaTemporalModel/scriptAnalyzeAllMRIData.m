
%% Housekeeping
clear
close all

%% Are we searching or not?
% Do we want to conduct a search for the fMRI data, or just use the p0
% values and make plots?
mriSearchFlag = true;

% Do we wish to use the monotonic constraint upon surround index in the
% search?
useMonotonicConstraint = false;

% Which set of parameters will we investigate in the bootstrap analysis?
modelType = 'bootV1';

% Which seed will we use to guide this bootstrap search
whichSeed = 'best';

% Verbose?
verbose = true;

% How many bootstrap resamplings of the data to conduct
nBoots = 1;

% Where we will save the temporal model results
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults');

%% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);

%% Load the Mt. Sinai data
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
postReceptoralPaths = {'midget','bistratified','parasol'};

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
mriTemporalModel.meta.postReceptoralPaths = postReceptoralPaths;
mriTemporalModel.meta.nBoots = nBoots;

paramCounts.unique = 5;
paramCounts.lgn = 1;
paramCounts.v1fixed = 1;
paramCounts.v1eccen = nEccs*2;
paramCounts.v1total = paramCounts.v1fixed+paramCounts.v1eccen;
mriTemporalModel.meta.paramCounts = paramCounts;

% Store a source version of the output variable.
mriTemporalModelSource = mriTemporalModel;

% Loop over subjects
for whichSub = [1 2]

    % Loop over bootstraps
    for bb = 1:nBoots

        % Get a resample with replacement of the acquisitions
        bootIdx = sort(datasample(1:nAcqs,nAcqs));

        % Assemble the data
        lgnY = []; lgnW = []; v1Y = []; v1W = [];
        for whichStim = 1:length(stimulusDirections)

            % Extract the relevant LGN data
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).lgn(bootIdx,:);
            lgnY = [lgnY mean(thisMatrix)];
            lgnW = [lgnW 1./std(thisMatrix)];

            % Extract the relevant V1 data acros eccentricities
            for ee = 1:nEccs
                thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['ecc' num2str(ee)])(bootIdx,:);
                v1Y = [v1Y mean(thisMatrix)];
                v1W = [v1W 1./std(thisMatrix)];
            end
        end

        % Load a search seed
        pMRI0 = storedSearchSeeds(whichSub,whichSeed);

        % Perform the search
        if mriSearchFlag

            % Report our progress
            curTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
            str=[curTime ' - Subject: ' subjects{whichSub} sprintf(', boot: %d ...',bb)];
            fprintf(str);

            % Let's see how long this takes
            tic

            % BADS it. We encounter occasional errors (parpool related?) so
            % place this in a try-catch.
            try
                result = fitMRIResponse(...
                    pMRI0,...
                    stimulusDirections,studiedEccentricites,studiedFreqs,...
                    v1Y,v1W,lgnY,lgnW,...
                    useMonotonicConstraint,modelType,verbose);
                result.pMRI0 = pMRI0;
            catch
                searchTimeSecs = toc();
                fprintf('error encountered. Skipping.\n');
                continue % skip this boot loop
            end

            % Report our search time and outcome
            searchTimeSecs = toc();
            str=[sprintf('fVal = %2.2f, search time (mins) = %2.1f',result.fVal,searchTimeSecs/60)  '\n'];
            fprintf(str);

        else
            result.pMRI = pMRI0;
            result.pMRI0 = pMRI0;
            result.fVal = nan;
        end

        % Report the parameters
        if verbose
            str = 'pMRI0 = [ ...\n';
            str = [str sprintf(repmat('%2.10f, ',1,paramCounts.unique),result.pMRI(1:paramCounts.unique)) ' ...\n'];
            str = [str sprintf(repmat('%2.10f, ',1,paramCounts.lgn*3),result.pMRI(paramCounts.unique+1:paramCounts.unique+paramCounts.lgn*3)) ' ...\n'];
            for ss=1:length(stimulusDirections)
                startIdx = paramCounts.unique+ paramCounts.lgn*3 + (ss-1)*paramCounts.v1total;
                str = [str sprintf(repmat('%2.10f, ',1,paramCounts.v1total),result.pMRI(startIdx+1:startIdx+paramCounts.v1total)) ' ...\n'];
            end
            str = [str ']; \n'];
            fprintf(str);
        end

        % Save the model parameters in an iteration structure and save this
        mriTemporalModel = mriTemporalModelSource;
        mriTemporalModel.(subjects{whichSub}).bootIdx = bootIdx;
        mriTemporalModel.(subjects{whichSub}).pMRI = result.pMRI;
        mriTemporalModel.(subjects{whichSub}).pMRI0 = result.pMRI0;
        mriTemporalModel.(subjects{whichSub}).fVal = result.fVal;
        mriTemporalModel.(subjects{whichSub}).v1Y = v1Y;
        mriTemporalModel.(subjects{whichSub}).v1W = v1W;
        mriTemporalModel.(subjects{whichSub}).lgnY = lgnY;
        mriTemporalModel.(subjects{whichSub}).lgnW = lgnW;

        % Save after each iteration, with a suffix that identifies the
        % particular resampling indicies
        bootLabel = regexprep(num2str(bootIdx),' +', '-');
        saveSpot = fullfile(saveDir,subjects{whichSub},['mriTemporalModel_' bootLabel '.mat']);
        save(saveSpot,'mriTemporalModel');

    end % boots

end % subjects



