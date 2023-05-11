
%% Housekeeping
clear
close all
verbose = true;
searchFlag = true;
nSearches = 4;
sortParamsPriorToSearch = false;
useMonotonicPenalty = true;

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
        if nn == 1
            p0 = [1.5 ...
                40.0000    1.6545    0.0470    8.0645    0.2500    1.1507   39.8275    1.7215    0.4402 ...
                40.0000    1.5158    0.0906    7.8748    0.3550    3.4444   37.5956    1.4989    0.6153 ...
                26.4518    1.1000    0.2256    7.2062    0.3622    5.9395   36.9056    1.3155    0.5311 ...
                12.1129    1.0000    0.2319    5.4367    0.3879    6.4357   34.5953    1.2660    0.5328 ...
                11.9229    0.9000    0.8838    5.2646    0.5467    6.0904    4.7764    1.1199    0.1937 ...
                11.4605    0.8000    3.6194    5.2897    0.5206    4.7633    2.8949    1.0106    1.5963 ...
                ];
        else
            p0 = results.(subjects{whichSub}).p;
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
        lb = [1.0 repmat([2.5, 0.25, 0.01],1,nCells*nEccs)];
        ub = [1.8 repmat([55, 2.5,  100],1,nCells*nEccs)];

        switch nn
            case 1
                % Only allow the gain and the Q value to vary
                lbG = p0; ubG = p0;
                idx = [1,4:nParams:length(lb)];
                lbG(idx) = lb(idx); ubG(idx) = ub(idx);
                lb = lbG; ub = ubG;
            case 2
                % Only allow the Q, gain, and filter frequency to vary
                lbG = p0; ubG = p0;
                idx = [1,4:nParams:length(lb),2:nParams:length(lb)];
                lbG(idx) = lb(idx); ubG(idx) = ub(idx);
                lb = lbG; ub = ubG;
            otherwise
                % use the full bounds
        end

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
        for ee=1:nEccs
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
            k = reshape(results.(subjects{whichSub}).p(2:end),nParams,nCells,nEccs);
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
% Add a penalty for non-monotonic ordering of corner frequencies and
% exponents across eccentricity. This was needed during the initial search
% through the parameter space

k=reshape(p(2:end),nParams,nCells,nEccs);

% For each parameter, calculate the penalty for two directions, and take
% the lower value

penalty = 0;
for cc = 1:3
    penalty = penalty + min([...
        sum(diff(squeeze(k(1,cc,:)))>0) ...
        sum(diff(squeeze(k(1,cc,:)))<0) ...
        ]);
    penalty = penalty + min([...
        sum(diff(squeeze(k(2,cc,:)))>0) ...
        sum(diff(squeeze(k(2,cc,:)))<0) ...
        ]);
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
    ttfAtEcc{ee} = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs,rgcTemporalModel);

end

% Reshape the responses into the dimension stim x ecc x freqs
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = ttfAtEcc{ee};
end

end


