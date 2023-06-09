% Here we fit the "sandwich" model to the data from each subject, and do so
% independently at each eccentricity.


%% Housekeeping
clear
close all
rng;
verbose = true;
sortParamsPriorToSearch = true;

% Define the source of the p0 values. By default, this is the result
% obtained by fitting the model separately at each eccentricity.
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
nSearches = 3;

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

% Turn off the warning
warning('off','bads:pbUnspecified');

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

    % Use the result obtained from fitting individual eccentricities
    p0 = results.(subjects{whichSub}).p;

    % Force the low-pass filter frequencies to be in descending order
    if sortParamsPriorToSearch
        k = reshape(p0(2:end),nParams,nCells,nEccs);
        for cc = 1:3
            k(1,cc,:) = sort(squeeze(k(1,cc,:)),'descend');
        end
        p0(2:end) = reshape(k,1,nParams*nCells*nEccs);
    end

    % Bounds on Q, corner frequency, exponentiation, gain. We put some
    % light upper bounds on the fits to prevent weird over-fitting that
    % happens for the low-amplitude chromatic responses at extreme
    % eccentricities.
    lb = [Q repmat([5 0.25 0.01 1 0.25 0.01 5 0.25 0.01],1,nEccs)];
    ub = [Q repmat([45 2.0 100 20 1.0 100 60 1.5 100],1,nEccs)];

    % Define the response function integrated across eccentricity
    myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs);
    myResponseMatrixInterp = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,interpFreqs);
    myResponseMatrixPlot = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,freqsForPlotting);
    myObj = @(p) norm(vectorize(Winterp).*(vectorize(Yinterp) - vectorize(myResponseMatrixInterp(p)))) + ...
        calcMonotonicPenalty(p,nParams,nCells,nEccs);

    % Loop over searches
    for nn = 1:nSearches

        % Start subsequent searches using the prior search result
        if nn>1
            p0 = p;
        end

        % BADS it
        [p,fVal] = bads(myObj,p0,lb,ub,[],[],[],optionsBADS);

        % Shut down the parpool to prevent memory leaks
        poolobj = gcp('nocreate');
        delete(poolobj);


    end % loop over searches

    % Get the response matrix and fVal
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
    save('cstResultsConstrained.mat','results')

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

% Restore the warnstate
warning(warnstate);



%% LOCAL FUNCTIONS

function penalty = calcMonotonicPenalty(p,nParams,nCells,nEccs)
% Add a penalty for non-monotonic ordering of corner frequencies and
% exponents across eccentricity.

k=reshape(p(2:end),nParams,nCells,nEccs);

penalty = 0;
for cc = 1:nCells
    penalty = penalty + min([ ...
        sum(diff(squeeze(k(1,cc,:)))>0) ...
        sum(diff(squeeze(k(1,cc,:)))<0) ...
        ]);
    penalty = penalty + min([ ...
        sum(diff(squeeze(k(2,cc,:)))>0) ...
        sum(diff(squeeze(k(2,cc,:)))<0) ...
        ]);
end

penalty = penalty / 1;

end


function response = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs)
% Assemble the response across eccentricity locations

% Fixed params of the analysis
nCells = 3;
nParams = 3;
blockLength = nParams*nCells;

% Loop over the passed eccentricities
parfor ee = 1:length(studiedEccentricites)

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


