
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
useMonotonicConstraint = true;


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
            pMRI0 = [ 0.9511927962, ...
                0.0121927595, 15.4630053043, 0.5059975982, 22.2882795334, 0.8433365583, 0.6106514931, 0.4058538914, 0.3140534639, 0.2288563490, 0.0193274736, 0.0198231936, 0.0210199356, 0.0407762527, 0.0458116531, 0.0555773973, 0.1293746233, ...
                0.0898599625, 25.2414739132, 0.3134130597, 9.9825036526, 0.5174748421, 0.5055159092, 0.4880894899, 0.3810717821, 0.1081254721, 0.0006752491, 0.4249371290, 0.4110758305, 0.3395687342, 0.2443841696, 0.1130095720, 0.1346800327, ...
                0.0492796206, 22.9782742262, 0.2441492319, 10.1444351673, 0.7291851044, 0.4346824884, 0.3814391375, 0.3215758085, 0.1540128469, 0.0622549772, 0.0796474218, 0.2485598326, 0.4388790131, 0.4008008242, 0.1855729818, 0.1770403385 ];
        case 2
            pMRI0 = [ 1.0090682446, ...
                0.0277694157, 12.3356672863, 0.4682257972, 20.8400594524, 0.5890172755, 0.5136457973, 0.3026875295, 0.0937087010, 0.0574547087, 0.0567557272, 0.0283756267, 0.0335035737, 0.0694713441, 0.0807542705, 0.2159285283, 0.9603156402, ...
                0.0874300769, 13.3532045800, 0.3060051915, 13.9739220623, 0.8312746300, 0.7818013132, 0.6472846015, 0.5454710624, 0.4042349860, 0.1931893200, 0.5614628771, 0.5980787728, 0.5849221377, 0.4465676009, 0.3343725225, 0.4180438078, ...
                0.0162586512, 14.9245809690, 0.1940235138, 13.5481184270, 0.2945510406, 0.2567878863, 0.2540839132, 0.0156124895, 0.0048501887, 0.0000000647, 0.0776062441, 0.2862449990, 0.6109366401, 0.5279773102, 0.6551479697, 0.9359545913 ];
    end

    % Perform the search
    if searchFlag

        % BADS it
        [pMRI,fVal] = fitMRIResponse(pMRI0,...
            stimulusDirections,studiedEccentricites,studiedFreqs,...
            v1Y,v1W,lgnY,lgnW,...
            useMonotonicConstraint);

    else
        pMRI = pMRI0;
        fVal = nan;
    end

    % Print the parameters in a format to be used as a seed in future searches
    str = ['pMRI0 = [ ' sprintf('%2.10f, ...\n',pMRI(1))];
    for ss=2:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-1,16)==0 && ss~=length(pMRI)
            str = [str '...\n'];
        end
    end
    str = [str(1:end-2) ' ];\n'];
    fprintf(str);

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
mriTemporalModel.meta.nFixedParams = 4;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 1;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
