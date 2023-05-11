
%% Housekeeping
clear
close all
verbose = false;
nSearches = 5;
sortParamsPriorToSearch = false;
useMonotonicPenalty = false;

Qvals = [0.5,0.75,1,1.25,1.5,1.75,2,2.25];

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

% A function to vectorize
vectorize = @(x) x(:);

% Turn off some warnings
warnstate = warning();

% Loop over Q vals
for qq = 1:length(Qvals)

    Q = Qvals(qq);

    % Loop over optimizations
    for nn = 1:nSearches

        % Loop over subjects
        for whichSub = 1:nSubs

            % Load the data
            Y = zeros(nSubs,nStims,nEccs,length(studiedFreqs));
            YstdSq = zeros(nSubs,nStims,nEccs,length(studiedFreqs));
            for eccIdx = 1:nEccs
                for whichStim = 1:nStims
                    thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(eccIdx)]);
                    Y(whichSub,whichStim,eccIdx,:) = mean(thisMatrix);
                    YstdSq(whichSub,whichStim,eccIdx,:) = std(thisMatrix).^2;
                end
            end
        end

        % Take the mean across subjects
        Y = squeeze(mean(Y,1));
        W = 1./sqrt(squeeze(sum(YstdSq,1))/nSubs);

        % We double the weights on the light-flux fits
        W(3,:,:) = W(3,:,:)*2;

        % If there was a prior analysis start with the last result
        if isfile('cstResults_minMod.mat')
            load('cstResults_minMod.mat','results')
        end

        % Use the prior solution as the p0 guess
        switch nn
            case 1
                p0 = [Q ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    120 1 0.1 120 1 10 120 1 1 ...
                    ];
            case 2
                p0 = [Q ...
                    120 1 0.1 20 1 10 20 1 1 ...
                    100 1 0.1 20 1 10 20 1 1 ...
                    80 1 0.1 20 1 10 20 1 1 ...
                    60 1 0.1 20 1 10 20 1 1 ...
                    40 1 0.1 20 1 10 20 1 1 ...
                    20 1 0.1 20 1 10 20 1 1 ...
                    ];
                p0(4:3:length(p0)) = results(qq).avgSub.p(4:3:length(p0));
            case 3
                p0 = [1.5 ...
  40.0000    1.6545    0.0470    8.0645    0.2500    1.1507   39.8275    1.7215    0.4402 ...
  40.0000    1.5158    0.0906    7.8748    0.3550    3.4444   37.5956    1.4989    0.6153 ...
   26.4518    1.1000    0.2256    7.2062    0.3622    5.9395   36.9056    1.3155    0.5311 ...
   12.1129    1.0000    0.2319    5.4367    0.3879    6.4357   34.5953    1.2660    0.5328 ...
   11.9229    0.9000    0.8838    5.2646    0.5467    6.0904    4.7764    1.1199    0.1937 ...
   11.4605    0.8000    3.6194    5.2897    0.5206    4.7633    2.8949    1.0106    1.5963 ...
                    ];
            otherwise
                p0 = results(qq).avgSub.p;
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
        switch nn
            case 1
                lb = [Q repmat([120, 1, 0.01],1,nCells*nEccs)];
                ub = [Q repmat([120, 1,  100],1,nCells*nEccs)];
            case 2
                lb = [Q repmat([2.5, 1, 0.01],1,nCells*nEccs)];
                ub = [Q repmat([120, 1,  100],1,nCells*nEccs)];
            otherwise
                lb = [Q repmat([2.5, 0.25, 0.01],1,nCells*nEccs)];
                ub = [Q repmat([120, 2.5,  100],1,nCells*nEccs)];
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
            myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites(ee),studiedFreqs,rgcTemporalModel);

            % Define the objective
            if nn>1 && useMonotonicPenalty
                myObj = @(p) norm(vectorize(W).*(vectorize(Y) - vectorize(myResponseMatrix(p)))) + ...
                    calcMonotonicPenalty(p,nParams,nCells,nEccs);
            else
                myObj = @(p) norm(vectorize(W(:,ee,:)).*(vectorize(Y(:,ee,:)) - vectorize(myResponseMatrix(p))));
            end

            % The place to define a non-linear constraint
            myNonbcon = [];

            % BADS it
            pSub = bads(myObj,p0(idx),lb(idx),ub(idx),[],[],myNonbcon,optionsBADS);

            % Store this eccentricity
            pPar{ee} = pSub;

        end

        % Sort out the loop p vals
        p=[];
        for ee = 1:nEccs
            idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];
            p(idx) = pPar{ee};
        end

        % Define the response function
        myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel);
        myObj = @(p) norm(vectorize(W).*(vectorize(Y) - vectorize(myResponseMatrix(p))));

        % Get the response matrix and fVal
        fVal = myObj(p);
        yFit = myResponseMatrix(p);

        % Report it
        fprintf('q = %2.2f, search = %d, fVal = %2.2f \n',Q,nn,fVal);

        % Store it
        results(qq).avgSub.p0 = p0;
        results(qq).avgSub.p = p;
        results(qq).avgSub.fVal = fVal;
        results(qq).avgSub.Y = Y;
        results(qq).avgSub.W = W;
        results(qq).avgSub.yFit = yFit;

        % Save it
        save('cstResults_minMod.mat','results')

        % Shut down the parpool to prevent memory leaks
        poolobj = gcp('nocreate');
        delete(poolobj);

    end

    % Plot it
    figure
    Y = results(qq).avgSub.Y;
    yFit = results(qq).avgSub.yFit;
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
    title(sprintf('q = %2.2f',Q))
    drawnow

    % Plot the params
    figure
    subLine = {'-',':'};
    paramNames = {'corner Freq','exponent','gain'};
    for pp = 1:nParams
        subplot(1,nParams,pp)
        for whichSub = 1:nSubs
            k = reshape(results(qq).avgSub.p(2:end),3,3,6);
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
    title(sprintf('q = %2.2f',Q))

end

% Restore the warnstate
warning(warnstate);



%% LOCAL FUNCTIONS

function penalty = calcMonotonicPenalty(p,nParams,nCells,nEccs)
% Add a penalty for non-monotonic ordering of corner frequencies and
% exponents across eccentricity. This was needed during the initial search
% through the parameter space

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

for ee = 1:length(studiedEccentricites)

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


