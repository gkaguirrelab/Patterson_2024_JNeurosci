function [response, responseMat, rfMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,studiedFreqs,rgcTemporalModel,paramCounts)

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
    activeCells = {}; lgnGain = []; rfPostRetinal = sym([]);

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Get the gain for this stimulus direction
        lgnGain = pMRI(paramCounts.unique + (ss-1)*paramCounts.lgn + 1);

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

            % Get the gain effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                modeledEccentricities(ee),stimulusContrastScale);

            % 2nd order low-pass fiter at the level of the retinto-
            % geniculate synapse
            rfPostRetinal = rfPostRetinal.*stageSecondOrderLP(lgnSecondOrderFc,lgnSecondOrderQ);

        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % Gain
        rfPostRetinal = (lgnGain/1e3)*rfPostRetinal;

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
