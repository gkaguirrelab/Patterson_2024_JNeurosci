
%% Housekeeping
clear
%close all


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
plotColor = {'r','b','k','g'};

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
        case 1 % fVal = 5.82753;
            pMRI0 = [...
                0.0104746437, 17.0244979858, 0.4969959259, 28.2221412659, 0.9507415771, 0.7042839050, 0.3849662781, 0.3405197144, 0.3355758667, 0.2324073792, 0.0127416611, 0.0145694733, 0.0337114334, 0.0407623291, 0.0580230713, 0.1623712540, ...
                0.1289903450, 16.8097305298, 0.4795207977, 29.3954753876, 0.6906654358, 0.2953659058, 0.1659790039, 0.0991832733, 0.0674797058, 0.0464511871, 0.1224637985, 0.5948043823, 0.8852195740, 0.6825754166, 0.2745433807, 0.2055805206, ...
                0.0682011032, 21.1593856812, 0.4451446533, 11.7552089691, 0.4788810730, 0.3627807617, 0.3428783417, 0.2408256531, 0.2327266693, 0.2152263641, 0.1135730743, 0.0561765671, 0.0885992050, 0.0381200790, 0.0991161346, 0.1839632034, ...
                0.0181715965, 21.8670654297, 0.4844245911, 11.7546081543, 0.8150615692, 0.5624813080, 0.5239078522, 0.1915817261, 0.1733913422, 0.1000473022, 0.0646680832, 0.0978424072, 0.0996274948, 0.1301408768, 0.0575757980, 0.0438117981 ...
                ];
        case 2 % fVal = 4.87029;
            pMRI0 = [...
                0.0255555247, 16.5456518834, 0.5028637315, 30.9526948805, 0.6457296962, 0.3783308995, 0.1422649911, 0.0133252568, 0.0130117555, 0.0033972142, 0.0146565018, 0.0196729188, 0.0439251990, 0.0630993931, 0.1767587359, 0.7513501852, ...
                0.0777374156, 16.5228143753, 0.5035740632, 31.7453588767, 0.8832296717, 0.0703800794, 0.0005334234, 0.0004254508, 0.0000745653, 0.0000697894, 0.0807823547, 0.6637335093, 0.9999961909, 0.9999993110, 0.8365181934, 0.8286270021, ...
                0.0750177029, 21.5815632860, 0.3423111681, 12.0955954917, 0.7320711422, 0.5172049241, 0.5151382418, 0.4824090677, 0.1885374246, 0.0210860621, 0.0585874533, 0.0733365953, 0.1283450486, 0.0784706920, 0.0688588814, 0.1144090461, ...
                0.0055124402, 21.5552491974, 0.3442547903, 11.9067317952, 0.8748793985, 0.8672860878, 0.4601877060, 0.1228676014, 0.1072251456, 0.0199814664, 0.0644168793, 0.1011525321, 0.0470758318, 0.1015484767, 0.1425107483, 0.3690340967 ...
                ];
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
    str = ['pMRI0 = [ ...\n'];
    for ss=1:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss,16)==0 && ss~=length(pMRI)
            str = [str '...\n'];
        end
    end
    str = [str(1:end-2) ' ... \n ];\n'];
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
mriTemporalModel.meta.nUniqueParams = 0;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
