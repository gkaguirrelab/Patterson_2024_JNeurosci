
%% Housekeeping
clear
close all
verbose = true;


%% Pick a subject, stimulus, and eccentricity band
whichSub = 1;


%% Load the data and basic model elements

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricitesFull = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};
nSubjects = 2;
nStims = 3;
nCells = 3;


% The number of acquisitions obtained for each measurement
nAcqs = 12;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

for whichSub = 1:nSubjects

    figure

    for eccIdx = 1:6
        nEccs = 1;
        studiedEccentricites = studiedEccentricitesFull(eccIdx);

        % Get the data
        for whichStim = 1:nStims
            for ee = length(eccIdx)
                whichEcc = eccIdx(ee);
                thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(whichEcc)]);
                Y(whichStim,ee,:) = mean(thisMatrix);
                W(whichStim,ee,:) = 1./std(thisMatrix);
            end
        end
        W(3,:,:) = W(3,:,:).*2;

        % Define the response function
        myResponseMatrix = @(p) returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel);

        % A function to vectorize
        vectorize = @(x) x(:);

        % Define the objective
        myObj = @(p) norm(vectorize(W).*(vectorize(Y) - vectorize(myResponseMatrix(p))));

        % No non-linear constraint at the moment
        myNonbcon = [];

        %% Set the bounds
        % Accelerating non-linearity
        % Q parameter
        % Vector of corner frequencies across eccentricity for each cell class
        % Vector of gains for each stimulus class
        lb =  [0.5,0.5,0.5,0.5,  ...
            repmat(5,1,nEccs), repmat(5,1,nEccs), repmat(10,1,nEccs), ...
            repmat(0.01,1,nCells*nEccs)];
        ub =  [2.5,2,2,2, ...
            repmat(30,1,nEccs), repmat(15,1,nEccs), repmat(45,1,nEccs), ...
            repmat(100,1,nCells*nEccs)];


        %% Set the p0
        switch subjects{whichSub}
            case 'gka'
                p0 =  [2,1,1,1, ...
                    repmat(30,1,nEccs), repmat(10,1,nEccs), repmat(35,1,nEccs), ...
                    repmat(0.1,1,nEccs), repmat(5,1,nEccs), repmat(1,1,nEccs)];
            case 'asb'
                p0 =  [1.5,1.5,1.5,1, ...
                    repmat(30,1,nEccs), repmat(10,1,nEccs), repmat(35,1,nEccs), ...
                    repmat(0.1,1,nEccs), repmat(5,1,nEccs), repmat(1,1,nEccs)];
        end

        % Options. Indicate that the objective function is deterministic, and
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

        % BADS it
        [p,fVal] = bads(myObj,p0,lb,ub,[],[],myNonbcon,optionsBADS)

        % Store it
        results.(subjects{whichSub}).p(eccIdx,:) = p;
        results.(subjects{whichSub}).fVal(eccIdx,:) = fVal;

        % Get the response
        yFit = myResponseMatrix(p);

        % Plot it
        for ss=1:nStims
            subplot(nStims,6,(ss-1)*6 + eccIdx)
            plot(squeeze(Y(ss,1,:)),'.k');
            hold on
            plot(squeeze(yFit(ss,1,:)),['-' plotColor{ss}]);
        end
        drawnow
    end
end
%     % Bode plot of the low-pass filter at the first eccentricity
% filterRF = stageFirstOrderLP(30,2);
% plotRF(filterRF);


foo = 1;


%% LOCAL FUNCTIONS

function [response, rfPostRetinal] = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel)


nEccs = length(studiedEccentricites);
nCells = 3;
nStims = length(stimulusDirections);

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

    % Loop over eccentricities
    for ee=1:nEccs

        ecc = studiedEccentricites(ee);

        % Loop over cell classes for this stimulus
        for cc=1:length(activeCells)

            % Get the scaling effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                ecc,stimulusContrastScale);

            % Extract the corner frequency parameter for this eccentricity
            % and cell class
            idx = 4 + (cellIdx(cc)-1)*nEccs + ee;
            cornerFrequency = p(idx);

            % Second order low pass filter at the level of retino-
            % geniculate synapse
            filter = stageSecondOrderLP(cornerFrequency,Q);
%            filter = filter ./ double(subs(filter,0));
            rfPostRetinal(cc) = rfPostRetinal(cc).*filter;

            % Obtain the exponent for the non-linear stage
            idx = 1+cellIdx(cc);
            cellExponent = p(idx);

            % Apply the non-linearity
            scaleVal = max(abs(double(subs(rfPostRetinal(cc),0:1:100))));
            rfPostRetinal = rfPostRetinal.^cellExponent;

            % Extract the gain for this cell class and stimulus direction
            idx = 4 + nCells*nEccs + (cellIdx(cc)-1)*nEccs + ee;
            gain = p(idx);

            % Gain
            rfPostRetinal(cc) = (gain/1000)*rfPostRetinal(cc).*(scaleVal/scaleVal^cellExponent);

        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % Derive the amplitude and phase from the Fourier model
        ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
        amplitude = abs(ttfComplex);
        phase = unwrap(angle(ttfComplex));

        if any(isnan(amplitude))
            foo=1;
        end

        response(ss,ee,:) = amplitude;
    end
end


end

