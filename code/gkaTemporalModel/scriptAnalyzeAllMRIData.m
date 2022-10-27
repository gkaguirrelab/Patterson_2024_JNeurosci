
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
        case 1 % fVal = 5.56342 
pMRI0 = [ 14.1315284371, 0.1637178987, 0.0144079423, 14.1365030408, 0.1213250220, 0.2083452475, 14.1370387375, 0.1237468243, 0.0154620746,  ...
19.5743161440, 0.6990626991, 29.5671971440, 0.9999556601, 0.9991493344, 0.4807035804, 0.4043660522, 0.1686772525, 0.0919928074, 0.8798409461, 0.8202997344, 1.9391235412, 2.0552867580, 3.2296952609, 4.6560890307, ...
19.5342242718, 0.5914926231, 29.5437857509, 0.8100783825, 0.1109152734, 0.1070137858, 0.0150517464, 0.0025852501, 0.0016614795, 0.7607503667, 2.9998657330, 4.1531111789, 2.8066742792, 1.0420949137, 0.8823529996, ...
19.5812153816, 0.2620808244, 10.9788304567, 0.6212356865, 0.3527286291, 0.2470120430, 0.2336025894, 0.0234890282, 0.0106283545, 3.0144434608, 0.8233963301, 2.9492289396, 4.3999650967, 3.7230858216, 6.2718994639, ...
41.5658599138, 0.4108411700, 10.9759245515, 0.7836344838, 0.5845842481, 0.5102277815, 0.4960044563, 0.1576134622, 0.1059410095, 4.3815265583, 5.9387723153, 4.9879750292, 3.6093874198, 1.8364786617, 1.1731457069 ... 
 ];
pMRI0 = [ 14.1509848833, 0.1524453074, 0.0153310165, 14.1652750969, 0.1528134942, 0.2051955241, 14.1552491486, 0.1528387696, 0.0152698299,  ...
19.4589212537, 0.6999720991, 29.7089674473, 0.9923491299, 0.9710890830, 0.4689561665, 0.4370403647, 0.3334224820, 0.1345343590, 0.9056850582, 0.8304520238, 1.9091819745, 2.0112036687, 2.5137119379, 5.6258597530, ...
19.4801053405, 0.6181894720, 29.7106980085, 0.5832179308, 0.1478988230, 0.1418548226, 0.0465576887, 0.0096996129, 0.0071535707, 0.7264291846, 2.8999377724, 4.1872875814, 2.8313885405, 1.0777787258, 0.8760118148, ...
19.4167846441, 0.2625908464, 10.7809215188, 0.3774427176, 0.3726565123, 0.2525582969, 0.2292559803, 0.0005589128, 0.0002837718, 2.2491661031, 0.8534269151, 2.7578709807, 4.4278748415, 4.0047659006, 6.5418000693, ...
38.7547850609, 0.4019058377, 10.7831016183, 0.8126747251, 0.5782607794, 0.5035249233, 0.4694375873, 0.2029336214, 0.0996335804, 4.7095033206, 6.0023219158, 5.1893657364, 3.7117403581, 1.5140301349, 0.6419357831 ... 
 ];

        case 2 % fVal = 4.5592
pMRI0 = [ 16.9521702826, 0.2564011373, 0.0289195460, 16.9151832163, 0.1098644547, 0.2028794373, 16.9332627207, 0.3312900424, 0.0138604907,  ...
14.3027577549, 0.4947960734, 28.7300854772, 0.9957950890, 0.6404848635, 0.2574268818, 0.0209027886, 0.0191291600, 0.0020987868, 0.6471415334, 0.8241390633, 2.0212030660, 2.7335790878, 7.7018873709, 32.7085615799, ...
14.3097254634, 0.2650850907, 28.7094716281, 0.4287352979, 0.1956231311, 0.0480252013, 0.0221532971, 0.0041389734, 0.0039771020, 2.3053495418, 6.0062291743, 10.0713253512, 8.5995535578, 6.1951971204, 5.8103362741, ...
14.3047544360, 0.3522058837, 10.4032430202, 0.5215743974, 0.2491518602, 0.2140900105, 0.1533369660, 0.0299713105, 0.0000824898, 0.3999196634, 0.5860188155, 0.3002529121, 1.7255136103, 5.8439547710, 17.1732794621, ...
44.6647769213, 0.2970009819, 10.4035698473, 0.9965207562, 0.8229722217, 0.4114404544, 0.0907435462, 0.0130804569, 0.0129990518, 3.5683417118, 4.6189411462, 5.5689502378, 4.9269433649, 3.2718934680, 1.9001256580 ... 
 ];
pMRI0 = [ 17.2246533632, 0.2789048910, 0.0296529210, 17.2647409141, 0.2789104342, 0.2092272544, 17.2050915658, 0.2793873489, 0.0148661762,  ...
14.1790354252, 0.4873219907, 28.6746282578, 0.9449516892, 0.7409727156, 0.2157115757, 0.0238312781, 0.0213572443, 0.0035259724, 0.6528011417, 0.8097884017, 1.9271248816, 2.7743697026, 7.5137959018, 31.0487744996, ...
14.1599076986, 0.2711931109, 28.6730350256, 0.4256471634, 0.1606073618, 0.0698900819, 0.0082490027, 0.0013934016, 0.0008855879, 1.7123265631, 5.9608129525, 10.4687383644, 9.0423364639, 7.4315749116, 5.4373528691, ...
14.1163358092, 0.3340698838, 10.0016410947, 0.6104362190, 0.2501534104, 0.1928754807, 0.1544524610, 0.0279316962, 0.0004897475, 0.7161051344, 0.5821384107, 0.3663670181, 2.2648384359, 6.0814075162, 15.8133978766, ...
42.0341461897, 0.2962514699, 10.0000001192, 0.9965842664, 0.9025683343, 0.4458149612, 0.0523285329, 0.0108442187, 0.0097243249, 3.9623665177, 4.5750525161, 5.7114482594, 4.3904520293, 2.9383908462, 1.3364374044 ... 
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
