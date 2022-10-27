
%% Housekeeping
clear
%close all


%% Create the RGC temporal sensitivity model
fitRGCFResponse


%% Are we searching or not?
% Do we want to conduct a search for the fMRI data, or just use the p0
% values and make plots?
searchFlag = false;

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
        case 1 % fVal = 5.52872
            pMRI0 = [ 14.1583633423, 0.1523025513, 0.0153379822, 14.1827392578, 0.1525039673, 0.2048143768, 14.1464233398, 0.1522041321, 0.0153530884,  ...
                19.5072937012, 0.6997261047, 29.7268524170, 0.9930160522, 0.9762893677, 0.4824035645, 0.4381652832, 0.4378402710, 0.2528305054, 0.9018505098, 0.8342977619, 1.9189728056, 2.0129641120, 2.5441992488, 6.6970722208, ...
                19.4789123535, 0.6118560791, 29.7302856445, 0.8240097046, 0.1464691162, 0.1408798218, 0.0457885742, 0.0089569092, 0.0077072144, 0.7242384275, 2.9272281747, 4.1935758857, 2.8448463195, 1.0329861914, 0.8764078753, ...
                19.4601440430, 0.2630630493, 10.7439880371, 0.3750350952, 0.3733139038, 0.2534576416, 0.2259094238, 0.0031341553, 0.0000030518, 2.2757646806, 0.8630999050, 2.8109164302, 4.4023788318, 3.9971332823, 6.5437913055, ...
                38.7123870850, 0.4010398865, 10.7332763672, 0.8226638794, 0.5757827759, 0.5027465820, 0.4667800903, 0.2063720703, 0.0990585327, 4.7103475355, 5.9975293882, 5.1802517214, 3.6984442502, 1.4331232506, 0.6401323995 ...
                ];

        case 2 % fVal = 4.62494
            pMRI0 = [ 17.2561692633, 0.2789886564, 0.0298463014, 17.2525762208, 0.2793695647, 0.2093741576, 17.2412980162, 0.2788229752, 0.0146849002,  ...
                14.1757942177, 0.4863741769, 28.6526598260, 0.9504602861, 0.7401523586, 0.2221389912, 0.0182465937, 0.0182063058, 0.0048223991, 0.6478636418, 0.8057045355, 1.9387292697, 2.7786515832, 7.5906455055, 31.4235298940, ...
                14.1719935276, 0.2687765557, 28.7128605992, 0.3265759848, 0.1611766674, 0.0715889439, 0.0093834564, 0.0004956920, 0.0000824004, 1.7505324276, 6.6208973856, 10.5466980312, 9.1119361998, 7.4528005102, 5.6280383488, ...
                14.1512143426, 0.3237100856, 10.0000203364, 0.6180180281, 0.5212171204, 0.2034285627, 0.1510728206, 0.0294813346, 0.0000383329, 0.7454038751, 1.0020953881, 0.3644301834, 2.2616777677, 6.0832686338, 15.8349247253, ...
                42.0911588334, 0.2956434639, 10.0004893541, 0.9978482295, 0.9010145575, 0.4567841277, 0.0606401913, 0.0100969899, 0.0084197845, 3.9579985144, 4.5549978290, 5.6847026516, 4.4416197277, 2.9160050042, 1.2039786845 ...
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
    str = sprintf('pMRI0 = [ %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f,  ...\n',pMRI(1:9));
    for ss=10:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-9,15)==0 && ss~=length(pMRI)
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
mriTemporalModel.meta.nFixedParams = 3;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 9;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
