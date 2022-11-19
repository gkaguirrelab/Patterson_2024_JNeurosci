function [response, responseMat, rfMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,paramCounts)

% Identify the studied eccentricities and stimulus frequencies
nCells = length(cellClasses);
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,nFreqs);

% Unpack the "unique" params
midgetChromSensitivityFactor = pMRI(1);
lgnSecondOrderFc = pMRI(2);
lgnSecondOrderQ = pMRI(3);
v1SecondOrderFc = pMRI(4);
v1SecondOrderQ = pMRI(5);

% Loop over eccentricities
parfor ee=1:nEccs

    activeCells = {};

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Unpack the V1 params (organized by stimulus)
        v1SurroundDelay = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + 1);
        v1SurroundIndex = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*0 + ee);
        v1Gain = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*1 + ee);

        % Identify which cell classes are relevant for this stimulus
        % direction
        switch stimulusDirections{ss}
            case 'LminusM'
                activeCells = {'midget'};
            case 'S'
                activeCells = {'bistratified'};
            case 'LMS'
                activeCells = {'midget','parasol'};
        end

        % The number of cell classes
        nCellsActive = length(activeCells);

        % Need to have an explicit loop length to keep parfor happy
        for cc=1:2

            if cc>nCellsActive
                continue
            end

            % Unpack the LGN params (organized by cell class)
            cellBlockIdx = find(strcmp(activeCells{cc},cellClasses));

            % Get the gain effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Apply the param adjustment to effect of L-< contrast upon
            % midgets. This covers the fact that we don't exactly know the
            % relative effectiveness of our stimulus contrast levels for
            % midget and parasol cells, and further will be realted to the
            % LM cone ratio in a way that I still need to spell out
            if strcmp(activeCells{cc},'midget') && strcmp(stimulusDirections{ss},'LMS')
                stimulusContrastScale = stimulusContrastScale * midgetChromSensitivityFactor;
            end

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                studiedEccentricites(ee),stimulusContrastScale);

            % 2nd order low-pass fiter at the level of the LGN
            rfPostRetinal = rfPostRetinal.*stageSecondOrderLP(lgnSecondOrderFc,lgnSecondOrderQ);

        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % Second order low pass filter at the level of V1
        rfPostRetinal = rfPostRetinal.*stageSecondOrderLP(v1SecondOrderFc,v1SecondOrderQ);

        % Delayed surround subtraction at V1
        rfPostRetinal = rfPostRetinal - v1SurroundIndex * rfPostRetinal.*stageDelay(v1SurroundDelay/1000);

        % Gain
        rfPostRetinal = (v1Gain/1000)*rfPostRetinal;

        % Store the RF
        rfMatrix(ss,ee) = rfPostRetinal;

        % Derive the amplitude and phase from the Fourier model
        ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
        v1Amplitude = abs(ttfComplex);
        v1Phase = unwrap(angle(ttfComplex));

        % Store this response
        responseMat(ss,ee,:) = v1Amplitude;

    end
end

response = [];
for ss = 1:nStims
    for ee=1:nEccs
        thisVec = squeeze(responseMat(ss,ee,:))';
        % Detect if we have a weird response and make it weirder to
        % penalize it.
        if thisVec(1)>thisVec(2)
            thisVec(1:3)=thisVec(1:3)*100;
        end
        response = [response thisVec];
    end
end

end