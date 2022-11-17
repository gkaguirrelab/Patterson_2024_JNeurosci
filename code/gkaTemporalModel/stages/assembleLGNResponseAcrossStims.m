function [response, responseMat, rfMatrix] = assembleLGNResponseAcrossStims(pMRI,cellClasses,stimulusDirections,studiedFreqs,rgcTemporalModel,paramCounts)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nFreqs = length(studiedFreqs);

% For the LGN model, we consider the average LGN response across a set of
% evenly spaced retinal eccentricities
modeledEccentricities = 1:5:81;
nEccs = length(modeledEccentricities);

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,nFreqs);

% Unpack the "unique" params
midgetStimSensitivityFactor = pMRI(1);

% Loop over eccentricities
parfor ee=1:nEccs

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Identify which cell classes are relevant for this stimulus direction,
        % as well as the parameters of the cortical, post-receptoral channel
        % second order filter
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

            % Unpack the LGN params (organized by cell class)
            cellBlockIdx = find(strcmp(activeCells{cc},cellClasses));

            lgnSurroundDelay = pMRI(paramCounts.unique + (cellBlockIdx-1)*paramCounts.lgn + 1);
            lgnSurroundIndex = pMRI(paramCounts.unique + (cellBlockIdx-1)*paramCounts.lgn + 2);
            lgnGain = pMRI(paramCounts.unique + (cellBlockIdx-1)*paramCounts.lgn + 3);

            % Get the gain effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirections{ss});

            % Apply the param adjustment to effect of LMS contrast upon
            % midgets. This covers the fact that we don't exactly know the
            % relative effectiveness of our stimulus contrast levels for
            % midget and parasol cells
            if strcmp(activeCells{cc},'midget') && strcmp(stimulusDirections{ss},'LMS')
                stimulusContrastScale = stimulusContrastScale * midgetStimSensitivityFactor;
            end

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                activeCells{cc},stimulusDirections{ss},rgcTemporalModel,...
                studiedEccentricites(ee),stimulusContrastScale);


            % Delayed surround subtraction at the LGN
            rfPostRetinal(cc) = rfPostRetinal(cc) - lgnSurroundIndex * rfPostRetinal(cc).*stageDelay(lgnSurroundDelay/1000);
        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % Gain
        rfPostRetinal = lgnGain*rfPostRetinal;

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
        thisVec = mean(squeeze(responseMat(ss,:,:)));
        response = [response thisVec];
    end
end

end
