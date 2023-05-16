% Here we fit the "sandwich" model to the data from each subject, and do so
% independently at each eccentricity.


%% Housekeeping
clear
close all
rng;
verbose = true;

% Define where we will save the results
resultFileName = 'cstResultsAvgROI.mat';

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);


% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

roiNames = {'lgn','v1_avg','v2v3_avg'};
nROIs = length(roiNames);

% Fixed features of the model
nCells = 3; nParams = 3;
nSearches = 0;

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
nFreqs = length(studiedFreqs);
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

% Turn off the warning
warning('off','bads:pbUnspecified');

% Loop over subjects
for whichSub = 1:nSubs

    for rr = 1:nROIs

        % Load the data for this subject
        Y = zeros(nStims,1,length(studiedFreqs));
        W = zeros(nStims,1,length(studiedFreqs));

        for whichStim = 1:nStims
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(roiNames{rr});
            Y(whichStim,1,:) = mean(thisMatrix);
            W(whichStim,1,:) = 1./std(thisMatrix);
        end

        % Fitting the LMS stimuli involves both the midget and parasol
        % populations; we give these weights a bit of extra oomph here to
        % encourage a good match to the data here
        W(3,:,:) = W(3,:,:)*1.5;

        % Interpolate the data and weights to the intermediate temporal
        % frequencies
        for whichStim = 1:nStims
            Yinterp(whichStim,1,:) = interp1(1:nFreqs,squeeze(Y(whichStim,1,:)),1:0.5:nFreqs);
            Winterp(whichStim,1,:) = interp1(1:nFreqs,squeeze(W(whichStim,1,:)),1:0.5:nFreqs);
        end

        % We set the p0 seed to have a declining corner frequency across
        % eccentricity, and a uniform value of the exponentiation term that
        % varies by cell class. The gain values are adjusted to be in the
        % appropriate range for each cell class.
        if isfile(resultFileName)
            p0 = results.(subjects{whichSub}).(roiNames{rr}).p;
            p = p0;
        else
            p0 = [Q 30 1.5 0.1 20 0.5 10 60 1 1 ];
            p = p0;
        end

        % Bounds on Q, corner frequency, exponentiation, gain. We put some
        % light upper bounds on the fits to prevent weird over-fitting that
        % happens for the low-amplitude chromatic responses at extreme
        % eccentricities.
        lb = [Q repmat([5 0.25 0.01 1 0.25 0.01 30 0.25 0.01],1,1)];
        ub = [Q repmat([30 2.0 100 30 2.0 100 60 2.0 100],1,1)];

        p0(p0<lb) = lb(p0<lb);
        p0(p0>ub) = lb(p0>ub);

        % Loop over searches
        for nn = 1:nSearches

            % Start subsequent searches using the prior search result
            if nn>1
                p0 = p;
            end

            % Define the response function
            myResponseMatrix = @(p) returnAvgResponse(p,stimulusDirections,interpFreqs,rgcTemporalModel);

            % Define the objective
            myObj = @(p) norm(vectorize(Winterp(:,1,:)).*(vectorize(Yinterp(:,1,:)) - vectorize(myResponseMatrix(p))));

            % BADS it
            p = bads(myObj,p0,lb,ub,[],[],[],optionsBADS);

        end % loop over searches

        % Define the response function integrated across eccentricity
        myResponseMatrix = @(p) returnAvgResponse(p,stimulusDirections,studiedFreqs,rgcTemporalModel);
        myResponseMatrixInterp = @(p) returnAvgResponse(p,stimulusDirections,interpFreqs,rgcTemporalModel);
        myResponseMatrixPlot = @(p) returnAvgResponse(p,stimulusDirections,freqsForPlotting,rgcTemporalModel);
        myObj = @(p) norm(vectorize(Winterp).*(vectorize(Yinterp) - vectorize(myResponseMatrixInterp(p))));

        % Get the response matrix and fVal
        fVal = myObj(p);
        yFit = myResponseMatrix(p);
        yPlot = myResponseMatrixPlot(p);

        % Store it
        results.(subjects{whichSub}).(roiNames{rr}).p0 = p0;
        results.(subjects{whichSub}).(roiNames{rr}).p = p;
        results.(subjects{whichSub}).(roiNames{rr}).fVal = fVal;
        results.(subjects{whichSub}).(roiNames{rr}).Y = Y;
        results.(subjects{whichSub}).(roiNames{rr}).W = W;
        results.(subjects{whichSub}).(roiNames{rr}).yFit = yFit;
        results.(subjects{whichSub}).(roiNames{rr}).yPlot = yPlot;

        % Save it
        save(resultFileName,'results')

        % Report it
        fprintf(['subject = ' subjects{whichSub} ', ROI = ' roiNames{rr} ' fVal = %2.2f \n'],fVal);

    end

    % Plot fits
    figure
    for rr = 1:nROIs
        Y = results.(subjects{whichSub}).(roiNames{rr}).Y;
        yPlot = results.(subjects{whichSub}).(roiNames{rr}).yPlot;
        for ss=1:nStims
            subplot(nStims,nROIs,(rr-1)*nStims+ss)
            plot(log10(studiedFreqs),squeeze(Y(ss,1,:)),'.k');
            hold on
            plot(log10(freqsForPlotting),squeeze(yPlot(ss,:)),['-' plotColor{ss}]);
            ylim([-2 6]);
            refline(0,0);
        end
        if rr == 1
            title(roiNames{rr})
        end
    end

    drawnow

end

% Plot params
figure
subLine = {'-','--'};
yLimSets = {[0 60],[0 2],[-1.25 2]};
paramNames = {'corner Freq','exponent','gain'};
for pp = 1:nParams
    subplot(1,nParams,pp)
    for whichSub = 1:nSubs
        vals = [];
        for ss = 1:nStims
            for rr = 1:nROIs
                k = reshape(results.(subjects{whichSub}).(roiNames{rr}).p(2:end),3,3,1);
                vals(rr) = squeeze(k(pp,ss,:));
            end
            if pp == 3
                plot(1:nROIs,log10(vals),[subLine{whichSub} 'o' plotColor{ss}]);
            else
                plot(1:nROIs,vals,[subLine{whichSub} 'o' plotColor{ss}]);
            end
            hold on
        end
    end
    title(paramNames{pp});
    ylim(yLimSets{pp});
    if pp == 2
        refline(0,1);
    end
    xlim([0.5 nROIs+0.5]);
    a = gca();
    a.XTickLabel = roiNames;
end

% Shut down the parpool to prevent memory leaks
% poolobj = gcp('nocreate');
% delete(poolobj);

% Restore the warnstate
warning(warnstate);



%% LOCAL FUNCTIONS



function response = returnAvgResponse(p,stimulusDirections,studiedFreqs,rgcTemporalModel)
% Assemble the response across eccentricity locations

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the V1 cortical bins
nEccs = 5;
studiedEccentricites = linspace(2,64,nEccs);

% Fixed params of the analysis
nCells = 3;
nParams = 3;

% Loop over the passed eccentricities
thisTTF = {};
parfor ee = 1:length(studiedEccentricites)

    % Obtain the response at this eccentricity
    thisTTF{ee} = returnTTFAtEcc(p,stimulusDirections,studiedEccentricites(ee),studiedFreqs,rgcTemporalModel);

end

% Obtain the sum across eccentricities
for ee = 1:length(studiedEccentricites)
    response(ee,:,:) = thisTTF{ee};
end
response = squeeze(mean(response,1));

% Detect if the response is not band pass and in that case make it a
% bad fit so that we avoid finding these solutions
for ss = 1:length(stimulusDirections)
    [~,idx] = max(response(ss,:));
    if idx < 4
        response(ss,:) = 100;
    end
end



end


