
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
        case 1 % fVal = 4.88295
pMRI0 = [ 18.5008511894, 0.4891848190, 27.0623925090, 16.9971477261, 0.4796920441, 29.1966505131, 22.4467143879, 0.1722127060, 12.0548238662, ...
0.0235109805, 0.9705113939, 0.6273803444, 0.3974060847, 0.2851224405, 0.3104442046, 0.5018055099, 0.0129172631, 0.0157150302, 0.0339501024, 0.0373058585, 0.0478791749, 0.1364880895, ...
0.1993744303, 0.7026468406, 0.2206110692, 0.1956990698, 0.0868973089, 0.0108395507, 0.0122635268, 0.1357933110, 0.6447451635, 0.9298858775, 0.6450356555, 0.2446142228, 0.1850322608, ...
0.0585864957, 0.2022314101, 0.1959049463, 0.3453165591, 0.3523901701, 0.2084339112, 0.1192959562, 0.0459756916, 0.0540069998, 0.1015824562, 0.0508040318, 0.0900819404, 0.1421393248, ...
0.1062309761, 0.8745291829, 0.8045935065, 0.7379196197, 0.2851119548, 0.0965183467, 0.3904407427, 0.1858331740, 0.1847644428, 0.1545300801, 0.2583987923, 0.0616510237, 0.0002126631 ... 
                ];
        case 2 % fVal = 4.96657
pMRI0 = [ 17.3533947580, 0.4830972036, 25.5861928686, 15.7240894996, 0.4124307889, 29.5230765641, 22.4467143879, 0.1722127060, 12.2921562847, ...
0.0539512290, 0.5779103007, 0.3897031777, 0.1657375749, 0.0010336112, 0.0057023626, 0.0398915298, 0.0182620657, 0.0222699828, 0.0483920336, 0.0628477848, 0.1760250759, 0.7757618385, ...
0.1273895105, 0.2595884629, 0.0005539540, 0.0005305786, 0.0000264674, 0.0000521105, 0.0111849029, 0.1081168371, 0.8051076760, 1.2603190776, 1.2091609456, 0.8997536928, 0.8526930709, ...
0.0917321921, 0.6821748169, 0.6722351824, 0.6130420759, 0.4775021104, 0.1652177572, 0.0073862565, 0.0803462392, 0.0847299215, 0.1455886648, 0.1243428409, 0.0981762623, 0.1291534651, ...
0.0008998580, 0.9830863469, 0.8226870943, 0.5120292419, 0.3166944749, 0.2153708917, 0.0022236756, 0.1507180059, 0.2324998414, 0.1883637634, 0.1895004336, 0.2180773662, 0.7691321721 ... 
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
    str = sprintf('pMRI0 = [ %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ...\n',pMRI(1:9));
    for ss=10:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-9,13)==0 && ss~=length(pMRI)
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
mriTemporalModel.meta.nFixedParams = 1;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = 9;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
