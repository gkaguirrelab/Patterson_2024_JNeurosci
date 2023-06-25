% Validate the fit for vertex 53477 for subject GKA


%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;
nBoots = 10;

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(studiedFreqs);
interpFreqs = logspace(log10(1),log10(100),501);
nAcqs = 12;

% Define some ROI sets
roiSet = {'LGN','V1','V2/V3','hV4','MT'};

%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the LGN ROI
tmpPath = fullfile(localDataDir,'retinoFiles','LGN_bilateral.dtseries.nii');
LGNROI = cifti_read(tmpPath); LGNROI = LGNROI.cdata;

% Load the MT ROI
tmpPath = fullfile(localDataDir,'retinoFiles','MT.dtseries.nii');
MTROI = cifti_read(tmpPath); MTROI = MTROI.cdata;

% fmincon Options. Indicate that the objective function is deterministic,
% and handle verbosity
options = optimoptions('fmincon');
options.Display = 'none';


% Set some bounds
LB = [0 1 0.5 0.5];
UB = [5 10 3 3];
p0A = [1.5 5 1.1 1.5];
p0B = [4 1.5 1.5 1];

% Anonymous functions for the search. The "modalPenalty" enforces that the
% interpolated response has a uni-modal distribution
myResp = @(p) watsonTemporalModel(p,studiedFreqs);
myRespInterp = @(p) watsonTemporalModel(p,interpFreqs);
modalPenalty = @(p) sum(sum(sign(diff(sign(diff(myRespInterp(p)))))) == 0)*1e3;

%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over bootstraps
    for bb = 1:nBoots

        % Get a sampling (with replacement) of the 12 acquisitions
        bootIdx = datasample(1:nAcqs,nAcqs);

        for rr = 1:length(roiSet)

            switch roiSet{rr}
                case 'LGN'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (LGNROI == 1)));
                case 'V1'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 1)));
                case 'V2/V3'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea > 1) .* (vArea < 4) ));
                case 'hV4'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 4) ));
                case 'MT'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (MTROI == 1)));
            end
            nGood(rr) = length(goodIdx);

            % Loop over the stimuli
            for whichStim = 1:nStims

                % Loop through the stimuli and frequencies and obtain the data
                rawVals = [];
                adjustedVals = [];

                for ff = 1:length(allFreqs)
                    subString = sprintf(['f%dHz_' stimulusDirections{whichStim}],allFreqs(ff));
                    idx = find(contains(stimLabels,subString));
                    rawVals(ff,:,:) = results.params(goodIdx,idx);
                end

                % Adjust the values for the zero frequency
                for ff = 2:length(allFreqs)
                    adjustedVals(ff-1,:,:) = rawVals(ff,:,:)-rawVals(1,:,:);
                end

                % Take the mean across voxels
                adjustedVals = squeeze(mean(adjustedVals,2));

                % Get the mean across (bootstrap resampled) acquisitions
                W = 1./std(adjustedVals(:,bootIdx),0,2)';
                Y = mean(adjustedVals(:,bootIdx),2)';

                % The weighted objective
                myObj = @(p) norm( W .* ( Y - myResp(p))) ...
                    + modalPenalty(p);

                % Fit it
                [pA, fValA] = fmincon(myObj,p0A,[],[],[],[],LB,UB,[],options);
                [pB, fValB] = fmincon(myObj,p0B,[],[],[],[],LB,UB,[],options);
                if fValA < fValB
                    p = pA;
                else
                    p = pB;
                end

                % Get the fit at the plotting frequencies
                yFit = myRespInterp(p);

                % Determine the peak frequency in the log domain
                peakFreq(whichStim,rr,bb) = log10(interpFreqs(yFit==max(yFit)));

            end % stimuli

        end % ROIs

    end % Bootstraps

    peakFreqSEM = 10.^std(peakFreq,0,3)
    peakFreqMean = 10.^mean(peakFreq,3)
    nGood


end
