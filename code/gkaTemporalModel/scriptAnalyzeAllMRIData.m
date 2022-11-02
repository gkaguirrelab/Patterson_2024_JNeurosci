
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

% How many bootstrap resamplings of the data to conduct
nBoots = 5;

% Where we will save the temporal model results
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults');


%% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);


%% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEcc = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k','g'};
postReceptoralPaths = {'midget.LminusM','bistratified.S','parasol.LMS','midget.LMS'};

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
mriTemporalModel.meta.nFixedParams = 1;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 9;
mriTemporalModel.meta.nBoots = nBoots;

% Store a source version of the output variable.
mriTemporalModelSource = mriTemporalModel;

% Loop over subjects
for whichSub = [2 1]

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
            for ee = 1:nEcc
                thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['ecc' num2str(ee)])(bootIdx,:);
                v1Y = [v1Y mean(thisMatrix)];
                v1W = [v1W 1./std(thisMatrix)];
            end
        end

        % Load a search seed
        pMRI0 = storedSearchSeeds(whichSub);

        % Perform the search
        if mriSearchFlag

            % Report our progress
            curTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
            str=[curTime ' - Subject: ' subjects{whichSub} sprintf(', boot: %d ...',bb)];
            fprintf(str);

            % Let's see how long this takes
            tic

            % BADS it
            [pMRI,fVal] = fitMRIResponse(pMRI0,...
                stimulusDirections,studiedEccentricites,studiedFreqs,...
                v1Y,v1W,lgnY,lgnW,...
                useMonotonicConstraint);

            searchTimeSecs = toc();

            str=[sprintf('fVal = %2.2f, search time (mins) = %2.1f',fVal,searchTimeSecs/60)  '\n'];
            fprintf(str);

        else            
            pMRI = pMRI0;
            fVal = nan;
        end

        % Save the model parameters in an iteration structure and save this
        mriTemporalModel = mriTemporalModelSource;
        mriTemporalModel.(subjects{whichSub}).bootIdx = bootIdx;
        mriTemporalModel.(subjects{whichSub}).pMRI = pMRI;
        mriTemporalModel.(subjects{whichSub}).pMRI0 = pMRI0;
        mriTemporalModel.(subjects{whichSub}).fVal = fVal;
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




%     Print the parameters in a format to be used as a seed in future searches
%{
    nEcc = 6;
    nUniqueParams = mriTemporalModel.meta.nUniqueParams;
    nParamsPerBlock = mriTemporalModel.meta.nFloatByEccParams * nEcc + mriTemporalModel.meta.nFixedParams;
    pathIndex = 1;
    str = ['pMRI0 = [ ...\n' sprintf([repmat('%2.10f, ',1,nUniqueParams) '... '],pMRI(1:nUniqueParams)) '%% lgn \n'];
    for ss=nUniqueParams+1:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-nUniqueParams,nParamsPerBlock)==0 && ss~=length(pMRI)
            str = [str '... %% V1 ' postReceptoralPaths{pathIndex} ' \n'];
            pathIndex = pathIndex+1;
        end
    end
    str = [str(1:end-2) ' ... %% V1 ' postReceptoralPaths{pathIndex} ' \n ]; \n'];
    fprintf(str); 
%}