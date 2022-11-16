function [response, responseMat, rfMatrix] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*3;
nSubtractions = 2; % Perform two levels of surround-delayed subtraction for V1

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,2,nFreqs);

% Loop over eccentricities
parfor ee=1:nEccs

    % Loop through the stimulus directions and assemble the response
    for ss = 1:nStims

        % Identify which cell classes are relevant for this stimulus direction,
        % as well as the parameters of the cortical, post-receptoral channel
        % second order filter
        switch stimulusDirections{ss}
            case 'LminusM'
                cellClasses = {'midget'};
                preV1ParamBlockIndex = 1;
            case 'S'
                cellClasses = {'bistratified'};
                preV1ParamBlockIndex = 2;
            case 'LMS'
                cellClasses = {'midget','parasol'};
                preV1ParamBlockIndex = [1 3];
        end

        % The number of cell classes that are relevant for this stimulus
        % direction
        nCellClasses = length(cellClasses);

        for cc=1:2

            if cc>nCellClasses
                continue
            end

            % Get the gain effect of stimulus contrast
            stimulusContrastScale = returnStimulusContrastScale(cellClasses{cc},stimulusDirections{ss});

            % Apply the param adjustment to effect of LMS contrast upon
            % midgets. This covers the fact that we don't exactly know the
            % relative effectiveness of our stimulus contrast levels for
            % midget and parasol cells
            if strcmp(cellClasses{cc},'midget') && strcmp(stimulusDirections{ss},'LMS')
                stimulusContrastScale = stimulusContrastScale * pMRI(end);
            end

            % Get the post-retinal temporal RF
            rfPostRetinal(cc) = returnPostRetinalRF(...
                cellClasses{cc},stimulusDirections{ss},rgcTemporalModel,...
                studiedEccentricites(ee),stimulusContrastScale);

            % Grab the LGN parameters
            lgnSurroundDelay = pMRI((preV1ParamBlockIndex(cc)-1)*nParamsPerCellBlock + 1);
            lgnSurroundIndex = pMRI((preV1ParamBlockIndex(cc)-1)*nParamsPerCellBlock + 2);

            % Delayed surround subtraction at the LGN
            rfPostRetinal(cc) = rfPostRetinal(cc) - lgnSurroundIndex * rfPostRetinal(cc).*stageDelay(lgnSurroundDelay/1000);
        end

        % Add the postRetinal RFs together
        rfPostRetinal = sum(rfPostRetinal);

        % Grab the V1 MRI parameters that are fixed across
        % eccentricity
%         secondOrderFc = pMRI((ss-1)*nParamsPerCellBlock + 4);
%         secondOrderQ = pMRI((ss-1)*nParamsPerCellBlock + 5);
%% Short circuit -- just use the first param block because these are now matched
        secondOrderFc = pMRI(4);
        secondOrderQ = pMRI(5);

        % Grab the V1 MRI parameters that vary with eccentricity

        %% Short circuit the delay param to just use the value from the first eccen bin
        % This is because we have forced these params to match anyway.
        v1SurroundDelay = pMRI((ss-1)*nParamsPerCellBlock + nFixedParams + nEccs*0 + 1);

        %% Short circuit the delay param to just use the value from the first cell block
        % This is because we have forced these params to match across cell blocks.
%        v1SurroundIndex = pMRI((ss-1)*nParamsPerCellBlock + nFixedParams + nEccs*1 + ee);
        v1SurroundIndex = pMRI(nFixedParams + nEccs*1 + ee);
        v1Gain = pMRI((ss-1)*nParamsPerCellBlock + nFixedParams + nEccs*2 + ee);

        % Delayed surround subtraction at V1
        rfPostRetinal = rfPostRetinal - v1SurroundIndex * rfPostRetinal.*stageDelay(v1SurroundDelay/1000);

        % Second order low pass filter
        rfPostRetinal = rfPostRetinal.*stageSecondOrderLP(secondOrderFc,secondOrderQ);

        % Gain
        rfPostRetinal = v1Gain*rfPostRetinal;

        % Store the RF
        rfMatrix(ss,ee) = rfPostRetinal;

        % Derive the amplitude and phase from the Fourier model
        ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
        v1Amplitude = abs(ttfComplex);
        v1Phase = unwrap(angle(ttfComplex));

        % Store this response
        responseMat(ss,ee,1,:) = v1Amplitude;

    end
end

response = [];
for ss = 1:nStims
    for ee=1:nEccs
        thisVec = sum(squeeze(responseMat(ss,ee,:,:)));
        % Detect if we have a weird response and make it weirder to
        % penalize it.
        if thisVec(1)>thisVec(2)
            thisVec(1:3)=thisVec(1:3)*100;
        end
        response = [response thisVec];
    end
end

end