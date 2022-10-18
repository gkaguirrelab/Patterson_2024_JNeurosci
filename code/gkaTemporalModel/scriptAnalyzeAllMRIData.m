
%% Housekeeping
clear
close all


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
eccDegVals = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimuli = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};
modelType = {'chromatic','chromatic','luminance'};

% The number of parameters in the model that are fixed across eccentricity.
% We need this information later when we make a plot of the surround index
% values across ecccentricity.
nFixed = 5;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Loop over stimulus directions and subjects
for whichStim = 1:3
    for whichSub = 1:2

        % Extract the relevant LGN data
        lgnFreqX = studiedFreqs;
        lgnY = squeeze(LGN_V1mri(whichSub,whichStim,1,:,2))';
        lgnW = 1./(squeeze(LGN_V1mri(whichSub,whichStim,1,:,3))'-squeeze(LGN_V1mri(whichSub,whichStim,1,:,1))');

        % Extract the relevant V1 data acros eccentricities
        v1Eccentricity = []; v1FreqX = []; v1Y = []; v1W = [];
        for ee = 1:nEcc
            v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
            v1FreqX = [v1FreqX studiedFreqs];
            v1Y = [v1Y squeeze(V1ecc_mri(whichSub,whichStim,ee,:,2))'];
            v1W = [v1W 1./(squeeze(V1ecc_mri(whichSub,whichStim,ee,:,3))'-squeeze(V1ecc_mri(whichSub,whichStim,ee,:,1))')];
        end

        % Stored p0 vectors
        switch whichStim
            case 1 % L-M chromatic
                switch whichSub
                    case 1 % fVal = 2.969
                        p0 = [ 1.0089677959, 0.0051972309, 17.3851397584, 0.7622797110, 30.9843077824, 0.7730219744, 0.5649684331, 0.3045017661, 0.2167170675, 0.1390385038, 0.0899452484, 0.0070561856, 0.0133253718, 0.0278128534, 0.0308904811, 0.0392816115, 0.1052877128 ];
                    case 2 % fVal = 2.738
                        p0 = [ 0.9914443493, 0.0134409952, 15.7920980453, 0.4662187099, 25.6041669846, 0.4843533516, 0.4111606598, 0.2233021736, 0.0211288452, 0.0175668716, 0.0157412529, 0.0150189400, 0.0287246704, 0.0578832626, 0.0699152946, 0.1890430450, 0.8185415268 ];

                end

            case 2 % S chromatic
                switch whichSub
                    case 1 % fVal = 1.972
                        p0 = [ 1.0002533518, 0.0034581543, 15.8262230880, 0.4096923549, 28.6208783539, 0.9539320752, 0.5253737290, 0.3478449759, 0.1708656253, 0.0693765355, 0.0690653905, 0.0017579700, 0.0097139167, 0.0285618296, 0.0421547971, 0.0534671235, 0.1456647850 ];
                    case 2 % fVal = 2.184
                        p0 = [ 0.9998304266, 0.0029785313, 12.6641071355, 0.3240441408, 20.5018852511, 0.3427464756, 0.3426924917, 0.2213753842, 0.0103705501, 0.0094405520, 0.0089717728, 0.0044831517, 0.0174973165, 0.0572883245, 0.0824097819, 0.2266700785, 0.8227456771 ];
                end

            case 3 % LMS luminance
                switch whichSub
                    case 1 % fVal = 3.151
                        p0 = [ 0.9827876066, 0.0452700749, 27.1020242213, 0.1470951200, 11.7800635150, 0.5624464617, 0.5624463381, 0.5077721188, 0.3753404966, 0.0897572991, 0.0000000000, 0.1512488763, 0.2119570867, 0.2405013722, 0.2269220409, 0.1418322324, 0.1951357739 ];
                    case 2 % fVal = 1.813
                        p0 = [ 1.3325173640, 0.0444406984, 16.2963053747, 0.1550703974, 15.2305855792, 0.9999999053, 0.7767361900, 0.6187224112, 0.4671422073, 0.3040515537, 0.2052299644, 0.1424478589, 0.2918556743, 0.3687851181, 0.3536844623, 0.3878452496, 0.7192984705 ];
                end
        end

        % Perform the search
        if searchFlag

            % BADS it
            [p,fVal] = fitMRIResponse(p0,v1FreqX,v1Eccentricity,v1Y,v1W, ...
                lgnFreqX,lgnY,lgnW, ...
                modelType{whichStim}, useMonotonicConstraint  );

            % Print the parameters in a format to be used as a seed in future searches
            str = 'p0 = [ ';
            for ss=1:length(p); str = [str sprintf('%2.10f, ',p(ss))]; end
            str = [str(1:end-2) ' ];\n'];
            fprintf(str);

        else
            p = p0;
            fVal = nan;
        end

        % Save the model parameters and data
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).p = p;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).fVal = p;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1FreqX = v1FreqX;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1Eccentricity = v1Eccentricity;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1Y = v1Y;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1W = v1W;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.lgnFreqX = lgnFreqX;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.lgnY = lgnY;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.lgnW = lgnW;
        mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).modelType = modelType{whichStim};

    end

end

%% Save the temporalModel
mriTemporalModel.meta.studiedFreqs = studiedFreqs;
mriTemporalModel.meta.eccDegVals = eccDegVals;
mriTemporalModel.meta.subjects = subjects;
mriTemporalModel.meta.stimuli = stimuli;
mriTemporalModel.meta.plotColor = plotColor;
mriTemporalModel.meta.nFixed = nFixed;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
