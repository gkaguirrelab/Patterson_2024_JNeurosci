
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
nFixed = 4;

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
                    case 1
                        p0 = [ 0.0051733944, 17.7443456650, 0.7639210075, 31.0002222657, 0.7747612000, 0.6127590656, 0.2972972929, 0.1360071421, 0.1359263182, 0.0927698433, 0.0068886280, 0.0128428638, 0.0290055275, 0.0304355621, 0.0388186872, 0.0996935070 ];
                    case 2
                        p0 = [ 0.0133828630, 15.7521562429, 0.4666754463, 25.3015081618, 0.4830160737, 0.4262431235, 0.2295079715, 0.0174463252, 0.0162542588, 0.0153346117, 0.0145549174, 0.0292946630, 0.0587017634, 0.0699338339, 0.1885285911, 0.7710133852 ];
                end

            case 2 % S chromatic
                switch whichSub
                    case 1
                        p0 = [ 0.0035079706, 15.7265627384, 0.4072433710, 28.7392640114, 0.9477780581, 0.5269501925, 0.3438204050, 0.1653482199, 0.1611490965, 0.0593631983, 0.0018090010, 0.0099489689, 0.0283274651, 0.0411759615, 0.0547645092, 0.1861535311 ];
                    case 2
                        p0 = [ 0.0030536129, 12.6634385064, 0.3230330031, 19.2353888974, 0.3384617426, 0.3383702502, 0.2204579670, 0.0074615318, 0.0058358349, 0.0053904641, 0.0045622010, 0.0174996126, 0.0576356798, 0.0819787085, 0.2284944840, 0.8201328889 ];
                end

            case 3 % LMS luminance
                switch whichSub
                    case 1
                        p0 = [ 0.0452219531, 26.8281599277, 0.1463285612, 11.8611085859, 0.5624449326, 0.5624449286, 0.5048329259, 0.3765357863, 0.0946555463, 0.0000000514, 0.1520011654, 0.2135339521, 0.2417147120, 0.2286185037, 0.1440009575, 0.1976674318 ];
                    case 2
                        p0 = [ 0.0441551732, 18.2843490106, 0.1358198817, 15.3268209426, 0.9999982010, 0.7692314450, 0.6248496951, 0.4658081513, 0.3066235899, 0.2012351101, 0.1455077090, 0.2923153721, 0.3706124250, 0.3581364105, 0.3933179319, 0.7217680452 ];
                end
        end

        % Perform the search
        if searchFlag

            % BADS it
            [p,fVal] = fitMRIResponse(p0,v1FreqX, v1Eccentricity, v1Y, v1W, ...
                lgnFreqX, lgnY, lgnW, ...
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
