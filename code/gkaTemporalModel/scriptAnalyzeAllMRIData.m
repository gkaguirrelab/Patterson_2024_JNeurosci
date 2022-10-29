
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
nBoots = 2;

% Where we will save the temporal model results
savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');


%% Create the RGC temporal sensitivity model
fitRGCFResponse


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
postReceptoralPaths = {'midget.LminusM','parasol.LMS','bistratified.S','midget.LMS'};

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
mriTemporalModel.meta.nFixedParams = 0;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 11;
mriTemporalModel.meta.nBoots = nBoots;

% Loop over subjects
for whichSub = 1:2

    % Loop over bootstraps
    for bb = 1:nBoots

        % Get a resample with replacement of the acquisitions
        bootIdx = datasample(1:nAcqs,nAcqs);

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
        pMRI0 = storedSearchSeeds(whichSub,useMonotonicConstraint);

        % Perform the search
        if mriSearchFlag

            % Report our progress
            str=['Subject: ' subjects{whichSub} sprintf(', boot: %d ...',bb)];
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
    end

%     Print the parameters in a format to be used as a seed in future searches
%     nUniqueParams = 11;
%     pathIndex = 1;
%     str = ['pMRI0 = [ ...\n' sprintf('%2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ... ',pMRI(1:5)) '%% lgn \n'];
%     str = [str sprintf('%2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ... ',pMRI(6:11)) '%% V1 chromatic, achromatic \n'];
%     for ss=nUniqueParams+1:length(pMRI)
%         str = [str sprintf('%2.10f, ',pMRI(ss))];
%         if mod(ss-nUniqueParams,12)==0 && ss~=length(pMRI)
%             str = [str '... %% V1 ' postReceptoralPaths{pathIndex} ' \n'];
%             pathIndex = pathIndex+1;
%         end
%     end
%     str = [str(1:end-2) ' ... %% V1 ' postReceptoralPaths{pathIndex} ' \n ]; \n'];
%     fprintf(str);

    % Save the model parameters and data
    mriTemporalModel.(subjects{whichSub}).pMRI(:,bb) = pMRI;
    mriTemporalModel.(subjects{whichSub}).fVal(:,bb) = fVal;
    mriTemporalModel.(subjects{whichSub}).data.v1Y(:,bb) = v1Y;
    mriTemporalModel.(subjects{whichSub}).data.v1W(:,bb) = v1W;
    mriTemporalModel.(subjects{whichSub}).data.lgnY(:,bb) = lgnY;
    mriTemporalModel.(subjects{whichSub}).data.lgnW(:,bb) = lgnW;

    % Save after each iteration, in case something breaks during the search
    save(savePath,'mriTemporalModel');

end
