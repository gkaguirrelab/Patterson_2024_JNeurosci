% Here we fit the "sandwich" model to the data from each subject, and do so
% independently at each eccentricity.


%% Housekeeping
clear
close all
rng;
verbose = false;

% Define where we will save the results
resultFileName = 'cstResults.mat';

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
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Fixed features of the model
nCells = 3; nParams = 3;
nSearches = 1;

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

% If there was a prior analysis start with the last result
if isfile(resultFileName)
    load(resultFileName,'results')
end

% Loop over subjects
for whichSub = 1:nSubs

    % Load the data for this subject
    Y = zeros(nStims,nEccs,length(studiedFreqs));
    W = zeros(nStims,nEccs,length(studiedFreqs));
    for eccIdx = 1:nEccs
        for whichStim = 1:nStims
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(eccIdx)]);
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
            Yinterp(whichStim,eccIdx,:) = interp1(1:nEccs,squeeze(Y(whichStim,eccIdx,:)),1:0.5:nEccs);
            Winterp(whichStim,eccIdx,:) = interp1(1:nEccs,squeeze(W(whichStim,eccIdx,:)),1:0.5:nEccs);
        end
    end

    % We set the p0 seed to have a declining corner frequency across
    % eccentricity, and a uniform value of the exponentiation term that
    % varies by cell class. The gain values are adjusted to be in the
    % appropriate range for each cell class.
    p0 = [Q ...
        30 1.5 0.1 20 0.5 10 60 1 1 ...
        26 1.5 0.1 18 0.5 10 50 1 1 ...
        22 1.5 0.1 16 0.5 10 40 1 1 ...
        18 1.5 0.1 14 0.5 10 30 1 1 ...
        14 1.5 0.1 12 0.5 10 20 1 1 ...
        10 1.5 0.1 10 0.5 10 10 1 1 ...
        ];

    % Bounds on Q, corner frequency, exponentiation, gain. We put some
    % light upper bounds on the fits to prevent weird over-fitting that
    % happens for the low-amplitude chromatic responses at extreme
    % eccentricities.
    lb = [Q repmat([5 0.25 0.01 1 0.25 0.01 5 0.25 0.01],1,nEccs)];
    ub = [Q repmat([45 2.0 100 20 1.0 100 60 1.5 100],1,nEccs)];

    % Loop over searches
    for nn = 1:nSearches

        % Start subsequent searches using the prior search result
        if nn>1
            p0 = p;
        end

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
            [pSub,fVal] = bads(myObj,p0(idx),lb(idx),ub(idx),[],[],[],optionsBADS);

            % Store this eccentricity
            pPar{ee} = pSub;

            % Report the fVal
            fprintf('ecc: %d, fVal = %2.2f \n',ee,fVal);

        end

        % Sort out the loop p vals
        p=[];
        for ee = 1:nEccs
            idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];
            p(idx) = pPar{ee};
        end

    end % loop over searches

    % Define the response function integrated across eccentricity
    myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs);
    myResponseMatrixInterp = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,interpFreqs);
    myResponseMatrixPlot = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,freqsForPlotting);
    myObj = @(p) norm(vectorize(Winterp).*(vectorize(Yinterp) - vectorize(myResponseMatrixInterp(p))));

    % Get the response matrix and fVal
    fVal = myObj(p);
    yFit = myResponseMatrix(p);
    yPlot = myResponseMatrixPlot(p);

    % Store it
    results.(subjects{whichSub}).p0 = p0;
    results.(subjects{whichSub}).p = p;
    results.(subjects{whichSub}).fVal = fVal;
    results.(subjects{whichSub}).Y = Y;
    results.(subjects{whichSub}).W = W;
    results.(subjects{whichSub}).yFit = yFit;

    % Save it
    save(resultFileName,'results')

    % Report it
    fprintf(['subject = ' subjects{whichSub} ', fVal = %2.2f \n'],fVal);

    % Plot fits
    figure
    Y = results.(subjects{whichSub}).Y;
    yFit = results.(subjects{whichSub}).yFit;
    for ee=1:6
        for ss=1:nStims
            subplot(nStims,nEccs,(ss-1)*nEccs + ee)
            plot(log10(studiedFreqs),squeeze(Y(ss,ee,:)),'.k');
            hold on
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
    for pp = 1:nParams
        subplot(1,nParams,pp)
        k = reshape(results.(subjects{whichSub}).p(2:end),3,3,6);
        for ss = 1:nStims
            plot(squeeze(k(pp,ss,:)),[subLine{whichSub} plotColor{ss}]);
            hold on
        end
        title(paramNames{pp});
        if pp == 2
            refline(0,1);
        end
        ylim(yLimSets{pp});
    end

    drawnow

end

% Shut down the parpool to prevent memory leaks
poolobj = gcp('nocreate');
delete(poolobj);

% Restore the warnstate
warning(warnstate);



%% LOCAL FUNCTIONS



function response = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs)
% Assemble the response across eccentricity locations

% Fixed params of the analysis
nCells = 3;
nParams = 3;
blockLength = nParams*nCells;

% Loop over the passed eccentricities
for ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    startIdx = (ee-1)*blockLength + 1 + 1;
    subP = [p(1) p(startIdx:startIdx+blockLength-1)];

    % Obtain the response at this eccentricity
    thisTTF = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs);

    % Detect if the response is not band pass and in that case make it a
    % bad fit so that we avoid finding these solutions
    for ss = 1:length(stimulusDirections)
        [~,idx] = max(thisTTF(ss,:));
        if idx < 4
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


