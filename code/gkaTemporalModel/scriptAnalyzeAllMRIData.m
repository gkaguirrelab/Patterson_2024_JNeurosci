
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
pMRI0 = [ 1.0397891998, ...
0.0115309620, 16.3018083572, 0.5158085823, 21.4839649200, 0.8157885551, 0.5678926468, 0.3731405258, 0.3109553337, 0.0774698257, 0.0525119781, 0.0200757980, 0.0217995644, 0.0432195663, 0.0440402031, 0.0509014130, 0.1226553917, ...
0.0888103819, 24.3418288231, 0.3654591084, 10.2755260468, 0.8293621063, 0.4864381981, 0.4783443451, 0.3677952766, 0.1069352150, 0.0000001907, 0.3157110214, 0.3419089317, 0.3005995750, 0.2145528793, 0.1034460068, 0.1245760918, ...
0.0477205276, 24.7860074043, 0.2188711643, 10.0629281998, 0.6446970940, 0.4633342743, 0.3875329018, 0.3388053894, 0.2208882332, 0.1121688843, 0.0773582458, 0.2707262039, 0.4487543106, 0.4175052643, 0.1818003654, 0.2274374962 ];
        case 2
pMRI0 = [ 0.9962771246, ...
0.0275532112, 12.4023504714, 0.4638420162, 20.6931730310, 0.6092341829, 0.5134320378, 0.2992227216, 0.0885561513, 0.0547090205, 0.0966085547, 0.0298811570, 0.0339005290, 0.0695820302, 0.0807837107, 0.2174377013, 0.9452668657, ...
0.0874345897, 13.3601331188, 0.3134470137, 13.9376569875, 0.8880759695, 0.7818749943, 0.6651186018, 0.5469190507, 0.3920963850, 0.1745152335, 0.1, 0.5835490970, 0.5976162743, 0.4318443120, 0.3204144899, 0.4150938802, ...
0.0196661160, 14.9489997280, 0.1960036997, 13.2678723666, 0.0019429828, 0.2683633326, 0.2594631566, 0.0151868575, 0.0032029303, 0.0015617853, 0.0014312179, 0.2886372029, 0.6106402007, 0.5341180102, 0.6552043105, 0.9684086721 ];
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
