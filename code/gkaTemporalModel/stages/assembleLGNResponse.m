function [response, responseMat, rfMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,studiedFreqs,rgcTemporalModel,paramCounts,modelType)

% Identify the stimuli and stimulus frequencies
nStims = length(stimulusDirections);
nFreqs = length(studiedFreqs);

% For the LGN model, we consider the average LGN response across a set of
% evenly spaced retinal eccentricities
modeledEccentricities = 1:5:81;
nEccs = length(modeledEccentricities);

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,nFreqs);

% Unpack the "unique" params
n = pMRI(1);
lgnSecondOrderFc = pMRI(2);
lgnSecondOrderQ = pMRI(3);

% Loop over eccentricities
parfor ee=1:nEccs

    % Clear some variables to keep parfor happy
    activeCells = {}; rfPostRetinal = sym([]);
    lgnGain = []; lgnSurroundDelay = []; lgnSurroundIndex = [];

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Identify which cell classes are relevant for this stimulus
        % direction.
        switch stimulusDirections{ss}
            case 'LminusM'
                activeCells = {'midget'};
            case 'S'
                activeCells = {'bistratified'};
            case 'LMS'
                activeCells = {'midget','parasol'};
        end

        % The number of cell classes that are relevant for this stimulus
        % direction
        nCellsActive = length(activeCells);

        % Need to have an explicit loop length to keep parfor happy
        for cc=1:2

            if cc>nCellsActive
                continue
            end

            % Get the LGN gain, which can be modeled on a per-cell or
            % per-stimulus basis
            switch modelType
                case 'stimulus'
                    lgnGain = pMRI(paramCounts.unique + (ss-1)*paramCounts.lgn + 3);
                case {'cell','mix'}
                    cellIdx = find(strcmp(activeCells{cc},cellClasses));
                    lgnSurroundDelay = pMRI(paramCounts.unique + (cellIdx-1)*paramCounts.lgn + 1);
                    lgnSurroundIndex = pMRI(paramCounts.unique + (cellIdx-1)*paramCounts.lgn + 2);
                    lgnGain = pMRI(paramCounts.unique + (cellIdx-1)*paramCounts.lgn + 1);
                otherwise
                    error('not a specified model type')
            end

            % Get the gain effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                modeledEccentricities(ee),stimulusContrastScale);

            % 2nd order low-pass fiter at the level of the retinto-
            % geniculate synapse
            rfPostRetinal(cc) = rfPostRetinal(cc).*stageSecondOrderLP(lgnSecondOrderFc,lgnSecondOrderQ);

            % If we are in the cell model, perform the delayed, surround
            % subtraction for each cell class, prior to combination
            switch modelType
                case 'stimulus'
                case 'cell'
                    rfPostRetinal(cc) = rfPostRetinal(cc) - lgnSurroundIndex * rfPostRetinal(cc).*stageDelay(lgnSurroundDelay/1000);
            end

            % Gain
            rfPostRetinal(cc) = (lgnGain/1e3)*rfPostRetinal(cc);

        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % If we are in the stimulus-refered model, perform delayed surround
        % subtraction upon the response for this stimulus after combining
        % cell classes
        switch modelType
            case {'stimulus','mix'}
                lgnSurroundDelay = pMRI(paramCounts.unique + (ss-1)*paramCounts.lgn + 1);
                lgnSurroundIndex = pMRI(paramCounts.unique + (ss-1)*paramCounts.lgn + 2);
                rfPostRetinal = rfPostRetinal - lgnSurroundIndex * rfPostRetinal.*stageDelay(lgnSurroundDelay/1000);
            case 'cell'
        end

        % Store the RF
        rfMatrix(ss,ee) = rfPostRetinal;

        % Derive the amplitude and phase from the Fourier model
        ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
        v1Amplitude = abs(ttfComplex);
        v1Phase = unwrap(angle(ttfComplex));

        % Apply the non-linear transformation of neural to BOLD response
        v1Amplitude = v1Amplitude.^n;

        % Store this response
        responseMat(ss,ee,:) = v1Amplitude;

    end

end

% Assemble the response vector to return
response = [];
for ss = 1:nStims
    thisVec = mean(squeeze(responseMat(ss,:,:)));
    response = [response thisVec];
end
responseMat = squeeze(mean(responseMat,2));

end
