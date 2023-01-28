
%% Housekeeping
clear
close all


%% Initialize the random seed so that we get different perms each run
rng('shuffle');


%% Are we searching or not?
% Do we want to conduct a search for the fMRI data, or just use the p0
% values and make plots?
mriSearchFlag = true;

% Do we wish to use the monotonic constraint upon surround index in the
% search?
useMonotonicConstraint = false;

% Which cortical region to fit
corticalRegion = 'v1';

% What model type do we want? By cell or by stimulus?
%{
    modelTypes = {'stimulus','cell'};
%}
modelTypes = {'cell','stimulus'};

% Which set of parameters will we investigate in the bootstrap analysis?
%{
paramSearch = 'gainOnly';
paramSearch = 'noSurround';
paramSearch = 'full';
paramSearch = 'cortex';
%}
paramSearch = 'gainOnly';

% How many bootstrap resamplings of the data to conduct
nBoots = 1;

% Verbose?
verbose = true;

% Where we will save the temporal model results
saveDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','v1');

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
mriTemporalModel.meta.useMonotonicConstraint = useMonotonicConstraint;

% Store a source version of the output variable.
mriTemporalModelSource = mriTemporalModel;

% Loop over bootstraps
for bb = 1:nBoots

    % Get a resample with replacement of the acquisitions
    if nBoots==1
        bootIdx = 1:12;
    else
        bootIdx = sort(datasample(1:nAcqs,nAcqs));
    end

    % Loop over subjects
    for whichSub = [1 2]

        % Assemble the data
        lgnY = []; lgnW = []; cortexY = []; cortexW = [];
        for whichStim = 1:length(stimulusDirections)

            % Extract the  LGN data
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).lgn(bootIdx,:);
            lgnY = [lgnY mean(thisMatrix)];
            lgnW = [lgnW 1./std(thisMatrix)];

            switch paramSearch
                case 'avgROIs'
                    % Extract the avg V1 response
                    thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).([corticalRegion '_avg'])(bootIdx,:);
                    cortexY = [cortexY mean(thisMatrix)];
                    cortexW = [cortexW 1./std(thisMatrix)];
                otherwise
                    % Extract the V1 response across eccentricities
                    for ee = 1:nEccs
                        thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).([corticalRegion '_ecc' num2str(ee)])(bootIdx,:);
                        cortexY = [cortexY mean(thisMatrix)];
                        cortexW = [cortexW 1./std(thisMatrix)];
                    end
            end
        end

        % Loop over models
        for mm=1:length(modelTypes)

            % Load a search seed
            pMRI0 = storedSearchSeeds(whichSub,modelTypes{mm},paramSearch);

            % Report our progress
            curTime = char(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss'));
            str=[curTime ' - Subject: ' subjects{whichSub} ', ' modelTypes{mm} sprintf(', boot: %d ...',bb)];
            fprintf(str);

            % Let's see how long this takes
            tic

            % BADS it.
            results = fitMRIResponse(...
                pMRI0,...
                stimulusDirections,studiedEccentricites,studiedFreqs,...
                cortexY,cortexW,lgnY,lgnW,...
                modelTypes{mm},useMonotonicConstraint,paramSearch,verbose);
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

            % Save the model parameters in an iteration structure and save this
            mriTemporalModel = mriTemporalModelSource;
            for fn = fieldnames(results)'
                mriTemporalModel.(subjects{whichSub}).(fn{1}) = results.(fn{1});
            end

            % Save after each iteration, with a suffix that identifies the
            % particular resampling indicies
            bootLabel = regexprep(num2str(bootIdx),' +', '-');
            saveSpot = fullfile(saveDir,modelTypes{mm},subjects{whichSub},['mriTemporalModel_' bootLabel '_' paramSearch '.mat']);
            save(saveSpot,'mriTemporalModel');

        end % model types

    end % subjects

end % boots



