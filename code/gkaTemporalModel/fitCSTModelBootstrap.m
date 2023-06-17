% This is the first stage of model fitting. The fMRI data is fit separately
% at each eccentricity by the [RGC][non-linear] model. Three parameters are
% used at each eccentricity and for each stimulus direction. The search is
% conducted in parallel across eccentricities. An outer loop creates a
% boot-strap resampling across acquisitions. This allows us to obtain the
% mean and distribution of the parameter fit values.



%% Housekeeping
clear
close all
rng('shuffle'); % Want a random order so we get different bootstraps
verbose = false;

% Define where we will save the results
resultFileName = 'cstResultsBootstrap.mat';

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the V1 cortical bins
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};
faceAlpha = 0.1; % Transparency of the shaded error region
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Fixed features of the model
nAcqs = 12; nCells = 3; nParams = 3;

% Control and plotting features
nBoots = 0;


% The "quality" parameter of the low-pass filter. In initial analyses we
% found that using a quality ("Q") parameter of unity for the low-pass
% filter stage produces a good result for both subjects, so we lock that in
% analysis here across subjects and eccentricities.
Q = 1.0;

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
studiedFreqs = [2 4 8 16 32 64];
interpFreqs = logspace(log10(2),log10(64),11);
freqsForPlotting = logspace(0,2,50);

% BADs Options. Indicate that the objective function is deterministic, and
% handle verbosity
optionsBADS.UncertaintyHandling = 0;
if verbose
    optionsBADS.Display = 'iter';
else
    optionsBADS.Display = 'off';
end

% The optimization toolbox is currently not available for Matlab running
% under Apple silicon. Detect this case and tell BADS so that it doesn't
% issue a warning
V = ver;
if ~any(strcmp({V.Name}, 'Optimization Toolbox'))
    optionsBADS.OptimToolbox = 0;
end

% A function to vectorize
vectorize = @(x) x(:);

% Save the warning state
warnstate = warning();

% If there was a prior analysis load it
if isfile(resultFileName)
    load(resultFileName,'results')
else
    results = [];
end

% We set the p0 seed to have a declining corner frequency across
% eccentricity, and a uniform value of the exponentiation term that
% varies by cell class. The gain values are adjusted to be in the
% appropriate range for each cell class.
p0 = [Q ...
    40 2.0 0.1 15 0.5 02 60 1.0 0.5 ...
    34 1.7 0.5 13 0.5 05 50 0.9 0.8 ...
    28 1.4 0.9 11 0.5 10 40 0.8 1.1 ...
    22 1.1 1.2 09 0.5 10 30 0.7 1.4 ...
    16 0.8 1.6 07 0.5 05 20 0.6 1.7 ...
    10 0.5 2.0 05 0.5 05 10 0.5 2.0 ...
    ];

% Loop over bootstraps
for bb = 1:nBoots

    % Open the parpool (and do so silentlty within an evalc)
    evalc('gcp()');

    % Get a sampling (with replacement) of the 12 acquisitions
    bootIdx = datasample(1:nAcqs,nAcqs);

    % Loop over subjects
    for whichSub = 1:nSubs

        % Report that we are starting the search
        fprintf(['boot %d, subject = ' subjects{whichSub} ' ... '],bb);
        tic

        % Load the data for this subject
        Y = zeros(nStims,nEccs,length(studiedFreqs));
        W = zeros(nStims,nEccs,length(studiedFreqs));
        for eccIdx = 1:nEccs
            for whichStim = 1:nStims

                % Extract the data
                thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(eccIdx)]);

                % Get this bootstrap resampling
                thisMatrix = thisMatrix(bootIdx,:);

                % Get the mean and weight
                Y(whichStim,eccIdx,:) = mean(thisMatrix);
                W(whichStim,eccIdx,:) = 1./std(thisMatrix);
            end
        end

        % Fitting the LMS stimuli involves both the midget and parasol
        % populations; we give these weights a bit of extra oomph here to
        % encourage a good match to the data here
        W(3,:,:) = W(3,:,:)*1.5;

        % Interpolate the data and weights to the intermediate temporal
        % frequencies
        for eccIdx = 1:nEccs
            for whichStim = 1:nStims
                Yinterp(whichStim,eccIdx,:) = interp1(1:nFreqs,squeeze(Y(whichStim,eccIdx,:)),1:0.5:nFreqs);
                Winterp(whichStim,eccIdx,:) = interp1(1:nFreqs,squeeze(W(whichStim,eccIdx,:)),1:0.5:nFreqs);
            end
        end

        % Bounds on Q, corner frequency, exponentiation, gain. We put some
        % light upper bounds on the fits to prevent weird over-fitting that
        % happens for the low-amplitude chromatic responses at extreme
        % eccentricities.
        lb = [Q repmat([5 0.25 0.01 1 0.25 0.01 5 0.25 0.01],1,nEccs)];
        ub = [Q repmat([45 2.0 100 20 1.5 100 60 1.5 100],1,nEccs)];

        % Loop over eccentricities
        parfor ee = 1:nEccs

            % Turn off the warning
            warning('off','bads:pbUnspecified');

            % Clear pSub
            pSub = [];

            % The parameter indices to work on
            idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];

            % Define the response function
            myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites(ee),interpFreqs);

            % Define the objective
            myObj = @(p) norm(vectorize(Winterp(:,ee,:)).*(vectorize(Yinterp(:,ee,:)) - vectorize(myResponseMatrix(p))));

            % BADS it
            pSub = bads(myObj,p0(idx),lb(idx),ub(idx),[],[],[],optionsBADS);

            % Store this eccentricity
            pPar{ee} = pSub;

        end

        % Sort out the loop p vals
        p=[];
        for ee = 1:nEccs
            idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];
            p(idx) = pPar{ee};
        end

        % Define the response function integrated across eccentricity
        myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs);
        myResponseMatrixInterp = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,interpFreqs);
        myResponseMatrixPlot = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,freqsForPlotting);
        myObj = @(p) norm(vectorize(Winterp).*(vectorize(Yinterp) - vectorize(myResponseMatrixInterp(p))));

        % Get the response matrix and fVal
        fVal = myObj(p);
        yFit = myResponseMatrix(p);
        yPlot = myResponseMatrixPlot(p);

        % Report the result
        timeMins = round(toc()/60);
        fprintf(' fVal = %2.2f, time [mins] = %d \n',fVal,timeMins);

        % Store it
        if ~isfield(results,subjects{whichSub})
            results.(subjects{whichSub}).p0 = [];
            results.(subjects{whichSub}).p = [];
            results.(subjects{whichSub}).fVal = [];
            results.(subjects{whichSub}).Y = {};
            results.(subjects{whichSub}).W = {};
            results.(subjects{whichSub}).yFit = {};
            results.(subjects{whichSub}).bootIdx = {};
        end

        results.(subjects{whichSub}).p0 = p0;
        results.(subjects{whichSub}).p(end+1,:) = p;
        results.(subjects{whichSub}).fVal(end+1) = fVal;
        results.(subjects{whichSub}).Y{end+1} = Y;
        results.(subjects{whichSub}).W{end+1} = W;
        results.(subjects{whichSub}).yFit{end+1} = yFit;
        results.(subjects{whichSub}).bootIdx{end+1} = bootIdx;

        % Save it
        save(resultFileName,'results')

    end

    % Shut down the parpool to prevent pseudo memory leaks that arise from
    % the persistence of symbolic toolbox elements in the workers
    evalc('delete(gcp(''nocreate''))');

    % Restore the warnstate
    warning(warnstate);

end % loop over bootstraps


% Finish up and plot; loop over subjects
for whichSub = 1:nSubs

    % Get the mean and the 67% CI of the parameters
    p = median(results.(subjects{whichSub}).p);
    pIQR = iqr(results.(subjects{whichSub}).p);
    pMat = sort(results.(subjects{whichSub}).p);
    pLow = p-pIQR/2;
    pHi = p+pIQR/2;
    Y = [];
    for ii = 1:length(results.(subjects{whichSub}).Y)
        Y(:,:,:,ii) = results.(subjects{whichSub}).Y{ii};
    end
    yIQR = iqr(Y,4);
    Y = median(Y,4);
    yLow = Y - yIQR/2;
    yHi = Y + yIQR/2;
    yPlot = returnResponse(p,stimulusDirections,studiedEccentricites,freqsForPlotting);

    % Plot fits
    figure
    for ee=1:6
        for ss=1:nStims
            subplot(nStims,nEccs,(ss-1)*nEccs + ee)

            thisVec = squeeze(Y(ss,ee,:))';
            thisLow = squeeze(yLow(ss,ee,:))';
            thisHi = squeeze(yHi(ss,ee,:))';

                        % Add a patch for the error
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ thisLow, fliplr(thisHi) ],...
                plotColor{ss},'EdgeColor','none','FaceColor',plotColor{ss},'FaceAlpha',faceAlpha);
            hold on

            plot(log10(studiedFreqs),squeeze(Y(ss,ee,:)),'.k');
            plot(log10(freqsForPlotting),squeeze(yPlot(ss,ee,:)),['-' plotColor{ss}]);
            ylim([-2 6]);
            refline(0,0);
        end
    end

    % Plot params
    figure
    subLine = {'-','-'};
    yLimSets = {[0 60],[0 2],[0 20]};
    paramNames = {'corner Freq','exponent','gain'};
    subP = reshape(p(2:end),nParams,nCells,nEccs);
    subPlow = reshape(pLow(2:end),nParams,nCells,nEccs);
    subPhi = reshape(pHi(2:end),nParams,nCells,nEccs);
    for pp = 1:nParams
        subplot(1,nParams,pp)

        for cc = 1:nCells

            thisVec = squeeze(subP(pp,cc,:))';
            thisLow = squeeze(subPlow(pp,cc,:))';
            thisHi = squeeze(subPhi(pp,cc,:))';

            % Add a patch for the error
            patch(...
                [log10(studiedEccentricites),fliplr(log10(studiedEccentricites))],...
                [ thisLow, fliplr(thisHi) ],...
                plotColor{cc},'EdgeColor','none','FaceColor',plotColor{cc},'FaceAlpha',faceAlpha);
            hold on
            plot(log10(studiedEccentricites),thisVec,[subLine{whichSub} plotColor{cc}]);
        end
        title(paramNames{pp});
        if pp == 2
            refline(0,1);
        end
        if pp == 3
            a = gca();
            a.YScale = 'log';
        end
        ylim(yLimSets{pp});
    end

end



%% LOCAL FUNCTIONS

function response = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs)
% Assemble the response across eccentricity locations

% Fixed params of the analysis
nCells = 3;
nParams = 3;

% Loop over the passed eccentricities
for ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];
    subP = p(idx);

    % Obtain the response at this eccentricity
    thisTTF = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs);

    % Detect if the response is not band pass and in that case make it a
    % bad fit so that we avoid finding these solutions
    for ss = 1:length(stimulusDirections)
        [~,idx] = max(thisTTF(ss,:));
        if idx < length(studiedFreqs)/3
            thisTTF(ss,:) = 100;
        end
    end

    % Store this loop result
    ttfAtEcc{ee} = thisTTF;

end

% Pull the fits out of the par pool cell array and place in a matrix
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = ttfAtEcc{ee};
end

end


