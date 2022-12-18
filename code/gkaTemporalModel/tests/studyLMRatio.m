% Housekeeping
clear

% Model and search type
modelType = 'stimulus';
paramSearch = 'full';

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
subjects = mriFullResultSet.meta.subjects;
nEccs = length(studiedEccentricites);

ratioVals = [1,1.25,1.5,2,4,8];

useMonotonicConstraint = false;
verbose = true;

% Loop over subjects and calculate the chromaticRG / luminance gain ratio
for whichSub = 1:length(subjects)

    % Assemble the data for the veridical dataset
    lgnY = []; lgnW = []; v1Y = []; v1W = [];
    for whichStim = 1:length(stimulusDirections)

        % Extract the  LGN data
        thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).lgn;
        lgnY = [lgnY mean(thisMatrix)];
        lgnW = [lgnW 1./std(thisMatrix)];

        % Extract the V1 response across eccentricities
        for ee = 1:nEccs
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(ee)]);
            v1Y = [v1Y mean(thisMatrix)];
            v1W = [v1W 1./std(thisMatrix)];
        end
    end

    % The average fit parameters across bootstraps
    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);

    for rr = 1:length(ratioVals)

        pMRI0 = pMRI;
        pMRI0(3) = ratioVals(rr);

        results = fitMRIResponse(...
            pMRI0,...
            stimulusDirections,studiedEccentricites,studiedFreqs,...
            v1Y,v1W,lgnY,lgnW,...
            'stimulus',useMonotonicConstraint,'full',verbose);

        % Store the results
        fVals(whichSub,rr) = results.fVal;
        pMRIResults(whichSub,rr,:) = results.pMRI;
    end
end
