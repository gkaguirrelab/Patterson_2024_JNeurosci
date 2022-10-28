
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
        case 1 % fVal = 5.99829
            pMRI0 = [ 12.7580675839, 0.2216328046, 0.0168856905, 0.2104816004, 0.0149722404, ...
                18.8684524664, 0.6999911702, 30.0686127458, 0.9963638382, 0.9435378242, 0.5068374404, 0.4507874305, 0.4434221172, 0.1819887379, 0.8946662106, 0.7654971599, 1.7762641918, 2.0648994509, 2.5511173598, 8.5562782891, ...
                18.8442010192, 0.6998521237, 30.0691210762, 0.7592912366, 0.1867649663, 0.1330747944, 0.0410067496, 0.0311164304, 0.0062370192, 0.5733913517, 2.7787164856, 4.0837011260, 2.7827739227, 1.0140216770, 0.7594296577, ...
                30.9480675477, 0.2532132396, 11.6006054845, 0.3423563966, 0.3181057562, 0.2247956649, 0.1964084852, 0.0828647843, 0.0001146987, 1.9080502795, 0.7720181254, 3.1886918603, 4.0145837119, 3.5883536397, 5.3830829850, ...
                30.9826550815, 0.2534385401, 11.6091857558, 0.8380121732, 0.6852821854, 0.6168135543, 0.5353163713, 0.2501092948, 0.1433036830, 6.7097979889, 8.3206215065, 5.8135870486, 4.0761158695, 1.4571386651, 0.8132767411 ...
                ];

        case 2 % fVal = 5.15352
            pMRI0 = [ 16.2089990362, 0.2485006001, 0.0350449744, 0.1820268166, 0.0136066973, ...
                14.7453904539, 0.3728790826, 30.6337734475, 0.9442481441, 0.8388733754, 0.3181913388, 0.0585735547, 0.0427402568, 0.0389790779, 0.5922679148, 0.7457920931, 1.9252944459, 2.6744853923, 7.0403982790, 31.9966349554, ...
                14.7183917285, 0.3727581348, 30.6726693658, 0.4478943266, 0.0759750013, 0.0435445487, 0.0235008070, 0.0180011918, 0.0169679177, 1.4691099341, 5.4616782348, 9.7857256367, 8.1786943332, 6.8461134837, 6.0710777710, ...
                32.9761217220, 0.2994634280, 10.0002098360, 0.5844829737, 0.5084927562, 0.3062582923, 0.1397208327, 0.0067924874, 0.0010427661, 0.4230190524, 1.7509563575, 0.3023717112, 2.9898315741, 4.3274536994, 12.2830631308, ...
                32.9761753624, 0.2996966663, 10.0001952772, 0.9922090812, 0.9494553246, 0.5073047666, 0.2042106104, 0.0696501230, 0.0138868256, 4.0138723342, 3.6535203727, 5.2865923017, 3.2100032718, 2.8511104081, 0.7784588511 ...
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
    nUniqueParams = 5;
    str = sprintf('pMRI0 = [ %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ...\n',pMRI(1:nUniqueParams));
    for ss=nUniqueParams+1:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-nUniqueParams,15)==0 && ss~=length(pMRI)
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
mriTemporalModel.meta.nUniqueParams = nUniqueParams;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
