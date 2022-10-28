
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
        case 1 % fVal = 5.90963
pMRI0 = [ 12.5464713317, 0.2116947352, 0.0167900426, 0.2105895916, 0.0154203754, ...
18.8577326294, 0.6999999999, 30.0610623239, 0.9915535146, 0.9476464373, 0.5018377482, 0.4565475496, 0.4448919144, 0.1802900279, 0.8859907238, 0.8309913006, 1.7725390712, 2.0493918273, 2.6055744065, 7.6484465218, ...
18.8575687213, 0.6999999999, 30.0639345422, 0.7430719873, 0.1908426744, 0.1218562596, 0.0330760624, 0.0314096121, 0.0043323757, 0.5680895895, 2.8409259587, 4.0901918776, 2.8091486840, 0.9888429670, 0.7459213014, ...
30.9427770879, 0.2532662169, 11.5559864910, 0.3316778087, 0.3241227247, 0.2293201071, 0.1949098121, 0.0742711783, 0.0005454157, 1.6065611973, 0.7106470615, 3.1962381980, 3.9998896621, 3.5914312281, 5.3828170884, ...
30.9558766941, 0.2532263007, 11.5558482232, 0.8494745031, 0.8297756940, 0.6187002586, 0.5512243262, 0.2524970833, 0.1437745741, 6.7790848191, 8.4623984049, 5.7556269059, 4.0518535577, 0.7177601826, 0.8452096698 ... 
 ];

        case 2 % fVal = 5.15223
pMRI0 = [ 16.2093003467, 0.2488001861, 0.0350617899, 0.1817352678, 0.0135776155, ...
14.7432534397, 0.3729892664, 30.6371337175, 0.9457955658, 0.8381944150, 0.3163870528, 0.0557254925, 0.0434300095, 0.0339942038, 0.5925804336, 0.7467029499, 1.9297409522, 2.6714913127, 7.0129004494, 31.9623597384, ...
14.7166985273, 0.3729356468, 30.6712754071, 0.4465603560, 0.0740879804, 0.0242479771, 0.0229250982, 0.0173201084, 0.0170098841, 1.4774250355, 5.4601488729, 9.7903958884, 8.2011225072, 6.8431216235, 6.0249467237, ...
32.9790417850, 0.2994748011, 10.0010076165, 0.5830877721, 0.5077618480, 0.3102626935, 0.1442157522, 0.0068452016, 0.0005129576, 0.3100553186, 1.7575214265, 0.9244686171, 2.9879124361, 4.3221554412, 12.2437430131, ...
32.9777492583, 0.2996853791, 10.0034207106, 0.9914230391, 0.9494037628, 0.5097826332, 0.2062790245, 0.0690515116, 0.0170384496, 4.0181709602, 3.7453289146, 5.2903415433, 3.2072681347, 2.8292254856, 0.4359396961 ... 
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
