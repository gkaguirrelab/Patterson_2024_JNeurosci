function [response, responseMat, rfMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,paramCounts,modelType,activeCellsLMS)

if nargin<9
    activeCellsLMS = {'midget','parasol'};
end

% Identify the studied eccentricities and stimulus frequencies
nCells = length(cellClasses);
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,nFreqs);

% The compressive non-linearity for neural-->BOLD response
n = pMRI(1);

% The corner frequency of the retino-geniculo synaptic filter
lgnSecondOrderFc = pMRI(2);
lgnSecondOrderQ = 0.45; % Fixed at 0.45

% The corner frequency of the geniculo-striate synaptic filter
v1SecondOrderFc = pMRI(2);
v1SecondOrderQ = 0.45; % Fixed at 0.45

% The LMRatio
rgcTemporalModel.LMRatio = pMRI(3);

% Loop over eccentricities
parfor ee=1:nEccs

    % Clear some variables to keep parfor happy
    activeCells = {}; rfPostRetinal = sym([]);
    v1Gain = []; v1SurroundDelay = []; v1SurroundIndex = [];

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Identify which cell classes are relevant for this stimulus
        % direction
        switch stimulusDirections{ss}
            case 'LminusM'
                activeCells = {'midget'};
            case 'S'
                activeCells = {'bistratified'};
            case 'LMS'
                activeCells = activeCellsLMS;
        end

        % The number of cell classes
        nCellsActive = length(activeCells);

        % Need to have an explicit loop length to keep parfor happy
        for cc=1:2

            % Part of the machinery to keep parfor happy
            if cc>nCellsActive
                continue
            end

            % Which cell class is relevant here?
            cellIdx = find(strcmp(activeCells{cc},cellClasses));

            % Get the V1 gain, which can be modeled on a per-cell or
            % per-stimulus basis. Also, if we are modeling on a cell level,
            % get the surround delay and index for this cell population.
            switch modelType
                case 'stimulus'
                    v1Gain = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*1 + ee);
                case 'cell'
                    v1SurroundDelay = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (cellIdx-1)*paramCounts.v1total + 1);
                    v1SurroundIndex = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (cellIdx-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*0 + ee);
                    v1Gain = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (cellIdx-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*1 + ee);
            end

            % Get the scaling effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                studiedEccentricites(ee),stimulusContrastScale);

            % Second order low pass filter at the level of retino-
            % geniculate synapse
            rfPostRetinal(cc) = rfPostRetinal(cc).*stageSecondOrderLP(lgnSecondOrderFc,lgnSecondOrderQ);

            % Second order low pass filter at the level of genciulo-
            % cortical synapse
            rfPostRetinal(cc) = rfPostRetinal(cc).*stageSecondOrderLP(v1SecondOrderFc,v1SecondOrderQ);

            % If we are in the cell model, perform the delayed, surround
            % subtraction for each cell class, prior to combination
            switch modelType
                case 'stimulus'
                case 'cell'
                    rfPostRetinal(cc) = rfPostRetinal(cc) - v1SurroundIndex * rfPostRetinal(cc).*stageDelay(v1SurroundDelay/1000);
            end

            % Gain
            rfPostRetinal(cc) = (v1Gain/1000)*rfPostRetinal(cc);

        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % If we are in the stimulus-refered model, perform delayed surround
        % subtraction upon the response for this stimulus after combining
        % cell classes
        switch modelType
            case 'stimulus'
                v1SurroundDelay = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + 1);
                v1SurroundIndex = pMRI(paramCounts.unique + paramCounts.lgn*nCells + (ss-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs*0 + ee);
                rfPostRetinal = rfPostRetinal - v1SurroundIndex * rfPostRetinal.*stageDelay(v1SurroundDelay/1000);
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
    for ee=1:nEccs
        thisVec = squeeze(responseMat(ss,ee,:))';
        % We can encounter a situation where the parameters imply a rising
        % gain at ever lower frequencies, but this is invisble to the fit
        % as we have a set of discrete frequencies. We try to detect this
        % case where gain is rising at low frequencies, and make this weird
        % response even weirder so as to penalize it.
        if thisVec(1)>thisVec(2)
            thisVec(1:3)=thisVec(1:3)*100;
        end
        response = [response thisVec];
    end
end

end