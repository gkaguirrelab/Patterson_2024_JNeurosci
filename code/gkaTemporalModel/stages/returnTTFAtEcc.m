function [ttfAtEcc, rfPostRetinalByStimulus] = returnTTFAtEcc(p,stimulusDirections,eccentricity,studiedFreqs,rgcTemporalModel)
% For a given visual field eccentricity and a given set of model
% parameters, return the temporal transfer function at a specified set of
% temporal frequencies

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

    clear rfPostRetinal;

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

        % Derive the amplitude and phase from the Fourier model
        ttfComplex = double(subs(rfPostRetinal(cc),studiedFreqs));
        amplitude(cc,:) = abs(ttfComplex);

    end

    % Assemble the responses to be returned
    ttfAtEcc(ss,:) = sum(amplitude,1);
    rfPostRetinalByStimulus{ss} = rfPostRetinal;

end


end
