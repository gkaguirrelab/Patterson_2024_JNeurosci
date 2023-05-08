
%% Housekeeping
clear
close all
verbose = true;
searchFlag = false;
nSearches = 1;
sortParamsPriorToSearch = false;
useMonotonicPenalty = false;

%% Load the data and basic model elements

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);

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

% The number of acquisitions obtained for each measurement. Might want this
% if we are going to do some boot-strapping
nAcqs = 12;

% Fixed features of the model
nCells = 3; nParams = 3;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

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

% Loop over optimizations
for nn = 1:nSearches

    % Loop over subjects
    for whichSub = 1:nSubs

        % Load the data
        Y = zeros(nStims,nEccs,length(studiedFreqs));
        W = zeros(nStims,nEccs,length(studiedFreqs));
        for eccIdx = 1:nEccs
            for whichStim = 1:nStims
                thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(eccIdx)]);
                Y(whichStim,eccIdx,:) = mean(thisMatrix);
                W(whichStim,eccIdx,:) = 1./std(thisMatrix);
            end
        end

        % We double the weights on the light-flux fits
        W(:,3,:) = W(:,3,:)*2;

        % If there was a prior analysis start with the last result
        load('cstResults.mat','results')

        % Use the prior solution as the p0 guess
        if isfield(results,subjects{whichSub})
            p0 = results.(subjects{whichSub}).p;
        else
            p0 = results.avgSub.p;
        end

        % Force the low-pass filter frequencies to be in descending order
        if sortParamsPriorToSearch
            k = reshape(p0(2:end),nParams,nCells,nEccs);
            for cc = 1:3
                k(1,cc,:) = sort(squeeze(k(1,cc,:)),'descend');
                k(2,cc,:) = sort(squeeze(k(2,cc,:)),'descend');
            end
            p0(2:end) = reshape(k,1,nParams*nCells*nEccs);
        end

        % Bounds on Q, corner frequency, exponentiation, gain
        lb = [1.0 repmat([ 5, 0.4, 0.01],1,nCells*nEccs)];
        ub = [1.8 repmat([45, 2.5,  100],1,nCells*nEccs)];

        % Define the response function
        myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel);

        % A function to vectorize
        vectorize = @(x) x(:);

        % Define the objective
        if useMonotonicPenalty
            myObj = @(p) norm(vectorize(W).*(vectorize(Y) - vectorize(myResponseMatrix(p)))) + ...
                calcMonotonicPenalty(p,nParams,nCells,nEccs);
        else
            myObj = @(p) norm(vectorize(W).*(vectorize(Y) - vectorize(myResponseMatrix(p))));
        end

        % The place to define a non-linear constraint
        myNonbcon = [];

        % BADS it
        if searchFlag
            [p,fVal] = bads(myObj,p0,lb,ub,[],[],myNonbcon,optionsBADS);
        else
            p=p0;
            fVal = myObj(p0);
        end

        % Get the response matrix
        yFit = myResponseMatrix(p);

        % Store it
        results.(subjects{whichSub}).p0 = p0;
        results.(subjects{whichSub}).p = p;
        results.(subjects{whichSub}).fVal = fVal;
        results.(subjects{whichSub}).Y = Y;
        results.(subjects{whichSub}).W = W;
        results.(subjects{whichSub}).yFit = yFit;

        % Save it
        if searchFlag
            save('cstResults.mat','results')
        end

        % Plot it
        figure
        Y = results.(subjects{whichSub}).Y;
        yFit = results.(subjects{whichSub}).yFit;
        for ee=1:6
            for ss=1:nStims
                subplot(nStims,nEccs,(ss-1)*nEccs + ee)
                plot(squeeze(Y(ss,ee,:)),'.k');
                hold on
                plot(squeeze(yFit(ss,ee,:)),['-' plotColor{ss}]);
                ylim([-2 6]);
                refline(0,0);
            end
        end
        drawnow

        % Bode plot of the low-pass filter at the first eccentricity
        % filterRF = stageFirstOrderLP(30,2);
        % plotRF(filterRF);

    end

    % Shut down the parpool to prevent memory leaks
    poolobj = gcp('nocreate');
    delete(poolobj);

    % Plot the params
    figure
    subLine = {'-',':'};
    paramNames = {'corner Freq','exponent','gain'};
    for pp = 1:nParams
        subplot(1,nParams,pp)
        for whichSub = 1:nSubs
            k = reshape(results.(subjects{whichSub}).p(2:end),3,3,6);
            for ss = 1:nStims
                plot(squeeze(k(pp,ss,:)),[subLine{whichSub} plotColor{ss}]);
                hold on
            end
        end
        title(paramNames{pp});
        if pp == 2
            refline(0,1);
        end
    end

end


foo = 1;


%% LOCAL FUNCTIONS

function penalty = calcMonotonicPenalty(p,nParams,nCells,nEccs)
% Add a penalty for non-monotonic ordering of corner frequencies across
% eccentricity

k=reshape(p(2:end),nParams,nCells,nEccs);

penalty = 0;
for cc = 1:3
    penalty = penalty + sum(diff(squeeze(k(1,cc,:)))>0);
    penalty = penalty + sum(diff(squeeze(k(2,cc,:)))>0);
end

penalty = penalty / 5;

end



function response = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel)
% Assemble the response across eccentricity locations

nCells = 3;
nParams = 3;
blockLength = nParams*nCells;

parfor ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    startIdx = (ee-1)*blockLength + 1 + 1;
    subP = [p(1) p(startIdx:startIdx+blockLength-1)];

    % Obtain the response at this eccentricity
    responseAtEcc{ee} = returnResponseAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs,rgcTemporalModel);

end

% Reshape the responses into the dimension stim x ecc x freqs
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = responseAtEcc{ee};
end

end


function responseAtEcc = returnResponseAtEcc(p,stimulusDirections,eccentricity,studiedFreqs,rgcTemporalModel)
% Obtain the modeled response at a given eccentricity

nCells = 3;
nStims = length(stimulusDirections);

% The params in p are organized as:
% - Q of 2nd-order low-pass filter
% - Frequency (Hz) of the corner frequency of the LP filter, x nCells
% - exponent of the non-linearity, x nCells
% - gain, x nCells

Q = p(1);

for ss=1:nStims
    switch stimulusDirections{ss}
        case 'LminusM'
            activeCells = {'midget'};
            cellIdx = 1;
        case 'S'
            activeCells = {'bistratified'};
            cellIdx = 2;
        case 'LMS'
            cellIdx = [1 3];
            activeCells = {'midget','parasol'};
    end

    % Loop over cell classes for this stimulus
    for cc=1:length(activeCells)

        % Get the scaling effect of stimulus contrast
        stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

        % Get the post-retinal temporal RF
        rfPostRetinal(cc) = returnPostRetinalRF(...
            activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
            eccentricity,stimulusContrastScale);

        % Extract the corner frequency parameter for this eccentricity
        % and cell class
        idx = 1 + (cellIdx(cc)-1)*nCells + 1;
        cornerFrequency = p(idx);

        % Second order low pass filter at the level of retino-
        % geniculate synapse
        filter = stageSecondOrderLP(cornerFrequency,Q);
        rfPostRetinal(cc) = rfPostRetinal(cc).*filter;

        % Obtain the exponent for the non-linear stage
        idx = 1 + (cellIdx(cc)-1)*nCells + 2;
        cellExponent = p(idx);

        % Apply the non-linearity
        scaleVal = max(abs(double(subs(rfPostRetinal(cc),0:1:100))));
        rfPostRetinal = rfPostRetinal.^cellExponent;

        % Extract the gain for this cell class and stimulus direction
        idx = 1 + (cellIdx(cc)-1)*nCells + 3;
        gain = p(idx);

        % Gain (plus backing out the scale change induced by the
        % exponentiation)
        rfPostRetinal(cc) = (gain/1000)*rfPostRetinal(cc).*(scaleVal/scaleVal^cellExponent);

    end

    % Add the postRetinal RFs together
    rfPostRetinal = sum(rfPostRetinal);

    % Derive the amplitude and phase from the Fourier model
    ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
    amplitude = abs(ttfComplex);
    phase = unwrap(angle(ttfComplex));

    responseAtEcc(ss,:) = amplitude;
end


end

