
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
            pMRI0 = [ 15.9288823756, 0.6083154351, 25.9726798639, 0.1516927078, 14.9026674218, 0.2505744883, ...
                0.0112027053, 23.7058968481, 0.8324431428, 0.5374327309, 0.3548559682, 0.2684237201, 0.2625181242, 0.2594704626, 0.0154302684, 0.0174757942, 0.0360333866, 0.0400823090, 0.0428825967, 0.1427896078, ...
                0.0920840123, 10.1694170421, 0.6176564738, 0.6009916501, 0.5999319907, 0.4991354119, 0.2092694358, 0.0013894393, 0.4267012434, 0.4092196599, 0.3500360483, 0.2496089209, 0.1143553689, 0.1037592194, ...
                0.0343509232, 23.1617868546, 0.6592568866, 0.4254767758, 0.3195293757, 0.3140736776, 0.0092653986, 0.0004335992, 0.0357879187, 0.1823894435, 0.3435748420, 0.3309946983, 0.1715137635, 0.1709332832 ];
        case 2
            pMRI0 = [ 11.2686634064, 0.4953022003, 20.1125431061, 0.1476696014, 10.5944919586, 0.1993751526,  ...
                0.0265002632, 20.9670639038, 0.6496463776, 0.5436134338, 0.3248447418, 0.1680938721, 0.0763168335, 0.0263640594, 0.0306491852, 0.0339584351, 0.0707569122, 0.0895252228, 0.2264118195, 1.0451469421, ...
                0.0907616997, 14.2565727234, 0.8696304321, 0.7629463196, 0.6952968597, 0.6235652924, 0.4706710815, 0.2612464905, 0.5726108551, 0.5820198059, 0.5274314880, 0.3972663879, 0.3189525604, 0.3830776215, ...
                0.0166347122, 13.2548713684, 0.3126876831, 0.3007892609, 0.2783817291, 0.0617530823, 0.0000007629, 0.0000000222, 0.0925846100, 0.3144321442, 0.6643066406, 0.5752792358, 0.6955280304, 1.1358489990 ];
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
    str = ['pMRI0 = [ ' sprintf('%2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ...\n',pMRI(1:6))];
    for ss=3:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-6,14)==0 && ss~=length(pMRI)
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
mriTemporalModel.meta.nFixedParams = 2;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 6;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
