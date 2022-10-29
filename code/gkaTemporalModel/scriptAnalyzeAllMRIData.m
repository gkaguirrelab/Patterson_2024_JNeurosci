
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
postReceptoralPaths = {'midget.LminusM','parasol.LMS','bistratified.S','midget.LMS'};

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
            % With monotonic constraint fVal = 5.82592
            pMRI0 = [ ...
                12.4648683038, 0.1974727350, 0.0164027326, 0.2106851201, 0.0152402088, ... % lgn
                18.9056068749, 0.6992170933, 29.5885303054, 30.8525628313, 0.2512669932, 11.5980746459, ... % V1 chromatic, achromatic
                0.9827368426, 0.9120300663, 0.4969400622, 0.4499106804, 0.4480630859, 0.3364435161, 0.9036672344, 0.8256234561, 1.8416708310, 1.9929656375, 2.7471324655, 7.6978269817, ... % V1 midget.LminusM
                0.7627563487, 0.1745144378, 0.1158846485, 0.0376191932, 0.0218891094, 0.0000378104, 0.6411567369, 2.8474927867, 4.0892546352, 2.7996272701, 0.9869840094, 0.6953022358, ... % V1 parasol.LMS
                0.7378422979, 0.3158069915, 0.2279572132, 0.1805807863, 0.0000541426, 0.0000329369, 1.6121338449, 0.7553756972, 3.2567305381, 4.0199905026, 3.6067176424, 5.5387705678, ... % V1 bistratified.S
                0.8669756048, 0.8264595550, 0.8202314845, 0.5742930332, 0.2312828546, 0.1297603187, 6.8545561469, 8.6045502884, 5.6830520717, 4.2375483123, 0.6030753915, 0.8088297225 ... % V1 midget.LMS
                ];

            % No monotonic constraint fVal = 5.76695
            pMRI0 = [ ...
                12.6704280259, 0.1936034098, 0.0166686890, 0.2078928001, 0.0152037158, ... % lgn
                18.9986599889, 0.6994403859, 29.3387743575, 31.1777582043, 0.2505297671, 11.5931272018, ... % V1 chromatic, achromatic
                0.9877896137, 0.9504461296, 0.4995610749, 0.4108317946, 0.5164756881, 0.7598901888, 0.8982917235, 0.8197531349, 1.8391523969, 1.9445522027, 2.3978482385, 6.6920877900, ... % V1 midget.LminusM
                0.7767462943, 0.1539912210, 0.0725080436, 0.0008021737, 0.0143844036, 0.0000963652, 0.6146171799, 2.8423936121, 4.0898950312, 2.7643015158, 0.9044463396, 0.7019457353, ... % V1 parasol.LMS
                0.5552172473, 0.3682746950, 0.2122287818, 0.1949797501, 0.0050239174, 0.0001724310, 1.4382524359, 0.6789144151, 3.2429015864, 3.9747648066, 3.6253579210, 5.3672825379, ... % V1 bistratified.S
                0.8920294924, 0.8326425004, 0.8208651459, 0.5901412179, 0.2304970280, 0.1221651037, 6.8963758750, 8.4373242918, 5.5376829649, 4.4770057977, 0.6601439654, 0.8007963881 ... % V1 midget.LMS
                ];
        case 2 % fVal = 5.06805
            pMRI0 = [ ...
                16.2844085693, 0.2524520874, 0.0352671814, 0.1845377350, 0.0137868500, ... % lgn
                14.0267181396, 0.3765632629, 29.2088317871, 32.5357055664, 0.2924240112, 10.1692199707, ... % V1 chromatic, achromatic
                0.9482376099, 0.8843109131, 0.3134979248, 0.0595214844, 0.0543579102, 0.0396499634, 0.5906873900, 0.7611677517, 1.9104963523, 2.6660494330, 7.1537434097, 31.4711469073, ... % V1 midget.LminusM
                0.4114685059, 0.1229171753, 0.0284622192, 0.0157546997, 0.0134979248, 0.0096115112, 1.6252690515, 5.4416227115, 9.6952727481, 8.3157612340, 6.8787920336, 6.1634068175, ... % V1 parasol.LMS
                0.5270263672, 0.4863906860, 0.3171875000, 0.1815063477, 0.0107833862, 0.0092132568, 0.8641012011, 1.8575443071, 0.6660816961, 3.0096155351, 4.1238079782, 12.2128747934, ... % V1 bistratified.S
                0.9886901855, 0.9245483398, 0.5364425659, 0.2519424438, 0.0756256104, 0.0146011353, 3.6232031479, 3.7093756342, 5.3217486042, 3.2047468904, 3.0168675730, 0.4079852891 ... % V1 midget.LMS
                ];
            % No monotonic constraint fVal = 4.96154
            pMRI0 = [ ...
                16.4064025879, 0.2356292725, 0.0350131226, 0.1849819946, 0.0140127563, ... % lgn
                13.9810180664, 0.3741668701, 28.1854248047, 32.5036621094, 0.2931793213, 10.2984619141, ... % V1 chromatic, achromatic
                0.9970581055, 0.9469787598, 0.3069702148, 0.0427124023, 0.0024017334, 0.2813171387, 0.5960661608, 0.7677481589, 1.9422044977, 2.6412490862, 7.0745093140, 30.6166245868, ... % V1 midget.LminusM
                0.0469238281, 0.0414855957, 0.0036560059, 0.0121704102, 0.0000122070, 0.0069702148, 1.6284985910, 5.4594324220, 9.6499037772, 8.2201902595, 6.8154021779, 6.5010525443, ... % V1 parasol.LMS
                0.5409423828, 0.4930419922, 0.2440368652, 0.1217224121, 0.0031982422, 0.0000122070, 0.9447277542, 2.0803727042, 0.6742047092, 2.9881213255, 4.0613319779, 12.0162330907, ... % V1 bistratified.S
                0.9645568848, 0.9974975586, 0.5505920410, 0.2654785156, 0.0647766113, 0.0000122070, 3.6506762851, 3.5974531834, 5.3159554469, 3.2450834066, 2.7828123390, 0.4173585463 ... % V1 midget.LMS
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
    nUniqueParams = 11;
    pathIndex = 1;
    str = ['pMRI0 = [ ...\n' sprintf('%2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ... ',pMRI(1:5)) '%% lgn \n'];
    str = [str sprintf('%2.10f, %2.10f, %2.10f, %2.10f, %2.10f, %2.10f, ... ',pMRI(6:11)) '%% V1 chromatic, achromatic \n'];
    for ss=nUniqueParams+1:length(pMRI)
        str = [str sprintf('%2.10f, ',pMRI(ss))];
        if mod(ss-nUniqueParams,12)==0 && ss~=length(pMRI)
            str = [str '... %% V1 ' postReceptoralPaths{pathIndex} ' \n'];
            pathIndex = pathIndex+1;
        end
    end
    str = [str(1:end-2) ' ... %% V1 ' postReceptoralPaths{pathIndex} ' \n ]; \n'];
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
mriTemporalModel.meta.nFixedParams = 0;
mriTemporalModel.meta.nFloatByEccParams = 2;
mriTemporalModel.meta.nUniqueParams = nUniqueParams;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
save(savePath,'mriTemporalModel');
