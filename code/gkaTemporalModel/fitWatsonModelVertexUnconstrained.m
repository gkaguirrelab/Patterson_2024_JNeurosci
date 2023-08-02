
%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Place to save the fit results
savePath = '~/Desktop/VSS 2023/';

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
analysisIDs = { '64a72e7ac71249914655c704','64a72ed7e54b4460ea4bf010' };
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


%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Define the localSaveDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;


% The initial guess for the Watson parameter search
    p0A = [1.5 5 1.1 1.5];


%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Initialize the fit results
    fitResults = [];
    fitResults.p = cell(nVert,1);
    fitResults.fVal = cell(nVert,1);
    fitResults.rSquared = cell(nVert,1);
    fitResults.Y = cell(nVert,1);
    fitResults.yFitInterp = cell(nVert,1);
    fitResults.peakFreq = cell(nVert,1);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Fit any voxel with an R2 fit above the threshold in V1
    goodIdx = find(logical( (results.R2 > r2Thresh)  .* (vArea == 1) ));
    nGood = length(goodIdx);

    % Loop over the valid voxels
    parfor gg = 1:nGood

        % The idx of the goodIdx set that we will process
        thisIdx = goodIdx(gg);

        % Loop through the stimuli and frequencies and
        % obtain the data
        rawVals = [];
        adjustedVals = [];
        for dd = 1:nStims
            for ff = 1:length(allFreqs)
                subString = sprintf(['f%dHz_' stimulusDirections{dd}],allFreqs(ff));
                idx = find(contains(stimLabels,subString));
                rawVals(dd,ff,:) = results.params(thisIdx,idx);
            end

            % Adjust the values for the zero frequency
            for ff = 2:length(allFreqs)
                adjustedVals(dd,ff-1,:) = rawVals(dd,ff,:)-rawVals(dd,1,:);
            end
        end

        % Set our data and weights for this vertex
        W = 1./std(adjustedVals,0,3,"omitmissing");
        Y = mean(adjustedVals,3,"omitmissing");

        % Define some variables to keep par happy       
        p = []; peakFreq = []; peakAmp = [];
        fVal = []; p = []; yFitInterp = []; rSquared = [];

        % Fit each stimulus with the Watson TTF
        for dd = 1:nStims

            [p(dd,:),fVal(dd),~,yFitInterp(dd,:)] = fitWatsonModel(Y(dd,:),W(dd,:),studiedFreqs,p0A,interpFreqs,false)

            % Determine the proportion variance explained
            maxfVal = norm( W(dd,:) .* ( Y(dd,:) ));
            rSquared(dd) = (maxfVal - fVal(dd)) / maxfVal;

            % Determine the peak frequency
            peakFreq(dd) = mean(interpFreqs(yFitInterp(dd,:)==max(yFitInterp(dd,:))));

            % Determine the peak ampitude
            peakAmp(dd) = p(dd,1);

        end

        % Place the par loop results into a results file
        parLoop_idx{gg} = thisIdx;
        parLoop_p{gg} = p;
        parLoop_fVal{gg} = fVal;
        parLoop_rSquared{gg} = rSquared;
        parLoop_Y{gg} = Y;
        parLoop_yFitInterp{gg} = yFitInterp;
        parLoop_peakFreq{gg} = peakFreq;
        parLoop_peakAmp{gg} = peakAmp;

    end

    % Move the par loop results into the fit structure
    idxSet = cell2mat(parLoop_idx);

    fitResults.p(idxSet) = parLoop_p;
    fitResults.fVal(idxSet) = parLoop_fVal;
    fitResults.rSquared(idxSet) = parLoop_rSquared;
    fitResults.Y(idxSet) = parLoop_Y;
    fitResults.yFitInterp(idxSet) = parLoop_yFitInterp;
    fitResults.peakFreq(idxSet) = parLoop_peakFreq;
    fitResults.peakAmp(idxSet) = parLoop_peakAmp;

    % Save the results file for this subject
    filePath = fullfile(savePath,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFitUnconstrained_results.mat']);
    save(filePath,'fitResults')


end


