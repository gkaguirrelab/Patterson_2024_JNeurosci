
%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% These variables define the subject names and stimulus directions.
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

%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Initialize the fit results
    fitResults = [];
    % fitResults.p = cell(nVert,1);
    % fitResults.fVal = cell(nVert,1);
    % fitResults.rSquared = cell(nVert,1);
    % fitResults.Y = cell(nVert,1);
    % fitResults.yFitInterp = cell(nVert,1);
    % fitResults.peakFreq = cell(nVert,1);
    % fitResults.peakAmp = cell(nVert,1);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Fit any voxel with an R2 fit above the threshold
    goodIdx = find(logical( (results.R2 > r2Thresh) ));
    nGood = length(goodIdx);

    % Loop over the valid voxels
    nGood = 3;
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

            [p(dd,:),fVal(dd),~,yFitInterp(dd,:)] = fitWatsonModel(Y(dd,:),W(dd,:),studiedFreqs)

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
    
    fitResults.idxSet = idxSet;
    fitResults.p = parLoop_p;
    fitResults.fVal = parLoop_fVal;
    fitResults.rSquared = parLoop_rSquared;
    fitResults.Y = parLoop_Y;
    fitResults.yFitInterp = parLoop_yFitInterp;
    fitResults.peakFreq = parLoop_peakFreq;
    fitResults.peakAmp = parLoop_peakAmp;

    % Save the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFit_resultsTEST.mat']);
    save(filePath,'fitResults')

end


