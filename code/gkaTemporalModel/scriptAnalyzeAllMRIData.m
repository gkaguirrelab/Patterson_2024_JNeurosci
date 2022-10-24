
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
        case 1 % fVal = 5.62022;
            pMRI0 = [...
                0.0096514514, 16.9536849769, 0.4747080471, 26.7965140999, 0.9494504243, 0.6191432365, 0.4177474823, 0.3531881250, 0.3417461254, 0.2530609437, 0.0143037875, 0.0169951651, 0.0353062475, 0.0398557569, 0.0591600474, 0.1748526851, ...
                0.0331277172, 16.0568186596, 0.4729856647, 27.6761134195, 0.5386474117, 0.2818790306, 0.2347393862, 0.1311477966, 0.0484268737, 0.0452431654, 0.0214998691, 0.0968586350, 0.1821976619, 0.1850298278, 0.1175224401, 0.1527603936, ...
                0.0846126476, 6.4253153367, 0.4471653647, 13.2133343493, 0.3763139432, 0.3760847502, 0.3551857584, 0.2299043211, 0.2298624531, 0.2291562536, 0.2276456545, 0.0933583570, 0.1639667250, 0.0195262319, 0.1098261684, 0.1875010255, ...
                0.0190960927, 22.7917484134, 0.4432269503, 13.2306410919, 0.7854896054, 0.5489842018, 0.5480254104, 0.1130288128, 0.0636353105, 0.0361660384, 0.0720787121, 0.0980291105, 0.0840466849, 0.1462666140, 0.0702575022, 0.0882488368 ...
                ];
        case 2
            pMRI0 = [...
                0.0253141287, 15.9147694229, 0.5087891015, 29.7061152597, 0.5802653073, 0.4130430978, 0.1593366459, 0.0185170320, 0.0178725225, 0.0173733213, 0.0164125700, 0.0199421476, 0.0454289773, 0.0624221357, 0.1828509109, 0.7331448772, ...
                0.0111104835, 10.2922708590, 0.5032965596, 29.7465669036, 0.8273481682, 0.1046735679, 0.0734414369, 0.0132576306, 0.0115967198, 0.0111313524, 0.0124488591, 0.1297043145, 0.2934839642, 0.3700417159, 0.4390176350, 0.6739828124, ...
                0.0904325822, 19.8985592897, 0.5176063669, 12.3619210009, 0.5692114182, 0.4803806015, 0.4765017052, 0.3924764121, 0.0727870340, 0.0037571036, 0.0903773207, 0.0551044979, 0.1410946365, 0.0590082187, 0.0620185154, 0.1240275678, ...
                0.0083895321, 20.0007404979, 0.5172311201, 12.3318046780, 0.9031572211, 0.8406957445, 0.3766229134, 0.1277947850, 0.0283913997, 0.0077974803, 0.0684830148, 0.1081490519, 0.0520393756, 0.1070084072, 0.0916081604, 0.1547213824 ...
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
mriTemporalModel.meta.nUniqueParams = 0;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
