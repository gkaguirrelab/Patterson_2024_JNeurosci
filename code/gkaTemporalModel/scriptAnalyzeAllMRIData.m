
%% Housekeeping
clear
close all


%% Create the RGC temporal sensitivity model
fitRGCFResponse


%% Are we searching or not?
% Do we want to conduct a search for the fMRI data, or just use the p0
% values and make plots?
searchFlag = true;

% Do we wish to use the monotonic constraint upon surround index in the
% search?
useMonotonicConstraint = false;


%% Load the Mt. Sinai data
% This should be the V1 and LGN area bold fMRI signal mean, and 95% CI. The
% matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x area
% (LGN 1, V1 2) x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_lgn_V1_BOLD.mat');
load(loadPath,'LGN_V1mri');

% Load the Mt. Sinai data
% This should be the V1 across eccentricity bold fMRI signal mean, and 95%
% CI. The matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x
% eccentricity x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_V1_ecc_BOLD.mat');
load(loadPath,'V1ecc_mri');

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEcc = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};

% The number of parameters in the model that are fixed across eccentricity.
% We need this information later when we make a plot of the surround index
% values across ecccentricity.
nFixedParams = 4;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Loop over subjects
for whichSub = 1:2
    lgnY = []; lgnW = [];
    v1Y = []; v1W = [];

    % Assemble the data
    for whichStim = 1:length(stimulusDirections)

        % Extract the relevant LGN data
        lgnY = [lgnY squeeze(LGN_V1mri(whichSub,whichStim,1,:,2))'];
        lgnW = [lgnW 1./(squeeze(LGN_V1mri(whichSub,whichStim,1,:,3))'-squeeze(LGN_V1mri(whichSub,whichStim,1,:,1))')];

        % Extract the relevant V1 data acros eccentricities
        for ee = 1:nEcc
            v1Y = [v1Y squeeze(V1ecc_mri(whichSub,whichStim,ee,:,2))'];
            v1W = [v1W 1./(squeeze(V1ecc_mri(whichSub,whichStim,ee,:,3))'-squeeze(V1ecc_mri(whichSub,whichStim,ee,:,1))')];
        end
    end

    % Stored p0 vectors
    switch whichSub
        case 1
            pMRI0 = [ 0.8863671795, ...
                0.0134860409, 18.1483947486, 0.5000765584, 21.2016258389, 0.7580165491, 0.5383087665, 0.3365183949, 0.2541487083, 0.0933484405, 0.0053778589, 0.0212477073, 0.0211576670, 0.0423804373, 0.0422141328, 0.0490258262, 0.0937200040, ...
                0.0857685912, 22.1618765593, 0.4247020379, 10.4420156777, 0.9734663874, 0.4701987937, 0.4378517017, 0.3489384219, 0.0927484572, 0.0070574537, 0.2510217130, 0.3128278330, 0.2691659108, 0.2009188607, 0.0960811749, 0.1269634813, ...
                0.0642368834, 38.5026419163, 0.2075465798, 5.0266880542, 0.5982929200, 0.5187450737, 0.3820989922, 0.4244025961, 0.1496773288, 0.2841300383, 0.1093606800, 0.3502666131, 0.4347161576, 0.5029014573, 0.1471649706, 0.2088705599 ];
        case 2
            pMRI0 = [ 1.0046916008, ...
                0.0287112510, 12.4640452862, 0.4762243032, 20.5744516850, 0.6210548878, 0.5393566608, 0.2984446526, 0.1115529060, 0.0505505323, 0.1056245804, 0.0295580626, 0.0334004164, 0.0670379400, 0.0802596807, 0.2271013260, 0.9598188400, ...
                0.0847976470, 12.9985713959, 0.3310333252, 13.9861059189, 0.8938520670, 0.7728738546, 0.6631578445, 0.5737783432, 0.4232458353, 0.1840734005, 0.5074352026, 0.5580109358, 0.5770705938, 0.4260205030, 0.3118672371, 0.4155191183, ...
                0.0097696412, 14.7498250008, 0.1956867456, 13.6855995655, 0.6553084850, 0.2811421156, 0.2825497389, 0.0013937950, 0.0016610861, 0.0037735224, 0.0510922670, 0.2876125574, 0.6287343502, 0.5228723288, 0.6475820541, 0.9739536047 ];
    end

    % Perform the search
    if searchFlag

        % BADS it
        [pMRI,fVal] = fitMRIResponse(pMRI0,...
            stimulusDirections,studiedEccentricites,studiedFreqs,...
            v1Y,v1W,lgnY,lgnW,...
            useMonotonicConstraint);

        % Print the parameters in a format to be used as a seed in future searches
        str = 'pMRI0 = [ ';
        for ss=1:length(pMRI); str = [str sprintf('%2.10f, ',pMRI(ss))]; end
        str = [str(1:end-2) ' ];\n'];
        fprintf(str);

    else
        pMRI = pMRI0;
        fVal = nan;
    end

    % Save the model parameters and data
    mriTemporalModel.(subjects{whichSub}).pMRI = pMRI;
    mriTemporalModel.(subjects{whichSub}).fVal = fVal;
    mriTemporalModel.(subjects{whichSub}).data.v1Y = v1Y;
    mriTemporalModel.(subjects{whichSub}).data.v1W = v1W;
    mriTemporalModel.(subjects{whichSub}).data.lgnY = lgnY;
    mriTemporalModel.(subjects{whichSub}).data.lgnW = lgnW;

end


%% Save the temporalModel
mriTemporalModel.meta.studiedFreqs = studiedFreqs;
mriTemporalModel.meta.studiedEccentricites = studiedEccentricites;
mriTemporalModel.meta.subjects = subjects;
mriTemporalModel.meta.stimulusDirections = stimulusDirections;
mriTemporalModel.meta.plotColor = plotColor;
mriTemporalModel.meta.nFixedParams = nFixedParams;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
