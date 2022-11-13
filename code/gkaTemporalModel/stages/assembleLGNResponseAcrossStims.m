function [response, responseMat] = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);

% How to divide up the pMRI veector
nParamsPerCellBlock = nFixedParams+nEccs*3;

% Loop through the stimulus directions and assemble the response
responseMat = zeros(nStims,2,nFreqs);

% Loop over the stims
for ss = 1:nStims

    % Identify which cell classes are relevant for this stimulus direction
    switch stimulusDirections{ss}
        case 'LminusM'
            cellClasses = {'midget'};
            paramBlockIndex = 1;
        case 'S'
            cellClasses = {'bistratified'};
            paramBlockIndex = 2;
        case 'LMS'
            cellClasses = {'parasol','midget'};
            paramBlockIndex = [3 1];
    end

    amplitudeResponseForCell = zeros(1,length(studiedFreqs));
    for cc = 1:length(cellClasses)

        % Get the gain effect of stimulus contrast
        stimulusContrastScale = returnStimulusContrastScale(cellClasses{cc},stimulusDirections{ss});

        % Apply the param adjustment to effect of LMS contrast upon
        % midgets. This covers the fact that we don't exactly know the
        % relative effectiveness of our stimulus contrast levels for
        % midget and parasol cells
        if strcmp(cellClasses{cc},'midget') && strcmp(stimulusDirections{ss},'LMS')
            stimulusContrastScale = stimulusContrastScale * pMRI(end);
        end

            % Grab the LGN parameters
        lgnSurroundDelay = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + 1);
        lgnSurroundIndex = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + 2);
        lgnGain = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + 3);

        % Get the TTF lgnAmplitude,lgnPhase
        amplitudeResponseForCell = amplitudeResponseForCell + ...
            returnlgnTTF(cellClasses{cc},stimulusDirections{ss},...
            rgcTemporalModel,stimulusContrastScale,...
            lgnSurroundDelay,lgnSurroundIndex,lgnGain,studiedFreqs);

        % Store this response
        responseMat(ss,cc,:) = amplitudeResponseForCell;

    end

end

response = [];
for ss = 1:nStims
    response = [response sum(squeeze(responseMat(ss,:,:)))];
end

end
