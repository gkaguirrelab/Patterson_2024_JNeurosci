
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
            pMRI0 = [ 15.9274065495, 0.6127423353, 25.5026811734, 0.1507032402, 15.0418103114, 0.2498860143, ...
                0.0113773800, 23.7346208841, 0.8302463949, 0.5453916177, 0.3401152462, 0.2616039351, 0.2593839511, 0.2231574506, 0.0156264380, 0.0174700320, 0.0351499543, 0.0373617038, 0.0487524569, 0.1461898759, ...
                0.0913746736, 10.0937556475, 0.6130318776, 0.6122673139, 0.6059100717, 0.5082148641, 0.2273885503, 0.0043636665, 0.4282640219, 0.4103404284, 0.3503185064, 0.2543968111, 0.1136030853, 0.1087179855, ...
                0.0300520833, 24.8813867941, 0.7510970801, 0.3865511104, 0.3608078092, 0.2819710121, 0.0171730921, 0.0018294439, 0.0365351513, 0.1782244295, 0.3349829912, 0.3303052038, 0.1723011658, 0.1982410327 ];
        case 2
            pMRI0 = [ 11.2682056427, 0.4890800118, 19.9337792397, 0.1453173399, 10.6194257736, 0.1993612051, ...
                0.0260793924, 20.9223008156, 0.6442087889, 0.5270112753, 0.3247524738, 0.1707157850, 0.0788476229, 0.0630104303, 0.0301693678, 0.0345637798, 0.0706949234, 0.0888468027, 0.2196089029, 0.8828469515, ...
                0.0897450459, 14.3934082985, 0.8976598740, 0.7726311684, 0.7227302074, 0.6227104187, 0.4718355894, 0.3247650623, 0.5481529236, 0.5780435801, 0.5541118383, 0.3864396811, 0.3042883873, 0.4408514500, ...
                0.0160578442, 11.3116735220, 0.3455569029, 0.2988475561, 0.2653712988, 0.0763711929, 0.0001698494, 0.0000366688, 0.0884656906, 0.3210479021, 0.6605762243, 0.5761026144, 0.6747804880, 1.0260866880 ];
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
    for ss=7:length(pMRI)
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
