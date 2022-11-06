function [response, responseMat] = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,freqs,rgcTemporalModel)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nFreqs = length(freqs);

% Loop through the stimulus directions and assemble the response
responseMat = zeros(nStims,2,nFreqs);

% Grab the lgn surround delay and index
surroundDelay = pMRI(1);
surroundIndex = pMRI(2);

% Loop over the stims
for ss = 1:nStims

    % Identify which cell classes are relevant for this stimulus direction
    switch stimulusDirections{ss}
        case 'LminusM'
            cellClasses = {'midget'};
        case 'S'
            cellClasses = {'bistratified'};
        case 'LMS'
            cellClasses = {'parasol','midget'};
    end

    amplitudeResponseForCell = zeros(1,length(freqs));
    for cc = 1:length(cellClasses)

        % Grab the LGN parameters, which are organized by RGC class
        switch cellClasses{cc}
            case 'midget'
                gain = pMRI(3);
            case 'bistratified'
                gain = pMRI(4);
            case 'parasol'
                gain = pMRI(5);
        end

        % Get the TTF lgnAmplitude,lgnPhase
        amplitudeResponseForCell = amplitudeResponseForCell + ...
            returnlgnTTF(cellClasses{cc},stimulusDirections{ss},...
            rgcTemporalModel,surroundDelay,surroundIndex,gain,freqs);

        % Store this response
        responseMat(ss,cc,:) = amplitudeResponseForCell;

    end


end

response = [];
for ss = 1:nStims
    response = [response sum(squeeze(responseMat(ss,:,:)))];
end

end
