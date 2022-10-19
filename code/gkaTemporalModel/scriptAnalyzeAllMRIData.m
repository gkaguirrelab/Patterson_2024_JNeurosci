
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
eccDegVals = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirection = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};

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
                    case 1 % fVal = 3.05672
                        p0 = [ 1.0365915596, 0.0026905593, 17.3618295789, 0.7657664865, 29.6612483263, 0.7719629049, 0.5324484706, 0.2823590875, 0.1830462933, 0.1764508307, 0.1291090548, 0.0107749701, 0.0119211376, 0.0266638994, 0.0316739678, 0.0390020609, 0.1110406220 ];
                    case 2 % fVal = 2.82654
                        p0 = [ 0.9620939866, 0.0068856616, 15.5466822535, 0.4683602266, 22.9088720679, 0.4899014831, 0.4327882633, 0.2212427914, 0.0365395889, 0.0186453119, 0.0138830900, 0.0220885724, 0.0262116492, 0.0551872253, 0.0694075748, 0.2092486992, 0.9071476012 ];


                end

            case 2 % S chromatic
                switch whichSub
                    case 1 % fVal = 1.75303
                        p0 = [ 1.7418847380, 0.0060342470, 13.1310718500, 0.4334087523, 28.8873984445, 0.9517605482, 0.3264243647, 0.2851549809, 0.1859624245, 0.0000000411, 0.0000000000, 0.0173579260, 0.1125805771, 0.2120167996, 0.2163356691, 0.1274167379, 0.1476888563 ];
                    case 2 % fVal = 1.93398
                        p0 = [ 0.5885573703, 0.0058107024, 10.7616140716, 0.3563334082, 9.4470165891, 0.2629912164, 0.2629888290, 0.2629875506, 0.0263283828, 0.0000445717, 0.0000000000, 0.0645086176, 0.2573866223, 0.5534312684, 0.4439267575, 0.5387061765, 0.7720702997 ];
                end

            case 3 % LMS luminance
                switch whichSub
                    case 1 % fVal = 3.06177
                        p0 = [ 1.5767331655, 0.0163077946, 28.7629187983, 0.1913331321, 11.4262334691, 0.5691350720, 0.5691349682, 0.4543243383, 0.3279572777, 0.0583557407, 0.0000003352, 0.2815376798, 0.2615394731, 0.2530629654, 0.2151177143, 0.1233465001, 0.1715851089 ];
                    case 2 % fVal = 2.12163
                        p0 = [ 0.3000000000, 0.0164496758, 17.4715094288, 0.1415145190, 15.4795711217, 1.0000000000, 0.7920254851, 0.6172449107, 0.4706962189, 0.3365233145, 0.2496413411, 0.3531578105, 0.4506660876, 0.5133960296, 0.4464027402, 0.4570651637, 0.8371949125 ];
                end
        end

        % Perform the search
        if searchFlag

            % BADS it
            [p,fVal] = fitMRIResponse(p0,v1FreqX,v1Eccentricity,v1Y,v1W, ...
                lgnFreqX,lgnY,lgnW, ...
                stimulusDirection{whichStim}, useMonotonicConstraint  );

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
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).p = p;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).fVal = p;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1FreqX = v1FreqX;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Eccentricity = v1Eccentricity;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Y = v1Y;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1W = v1W;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnFreqX = lgnFreqX;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnY = lgnY;
        mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnW = lgnW;

    end

end

%% Save the temporalModel
mriTemporalModel.meta.studiedFreqs = studiedFreqs;
mriTemporalModel.meta.eccDegVals = eccDegVals;
mriTemporalModel.meta.subjects = subjects;
mriTemporalModel.meta.stimulusDirection = stimulusDirection;
mriTemporalModel.meta.plotColor = plotColor;
mriTemporalModel.meta.nFixed = nFixed;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
