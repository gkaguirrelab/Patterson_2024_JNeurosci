function [response, responseMat] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*3;
nSubtractions = 2; % Perform two levels of surround-delayed subtraction for V1

% Initialize the response matrix
responseMat = zeros(nStims,nEccs,2,nFreqs);

% Loop through the stimulus directions and assemble the response
for ss = 1:nStims

    % Identify which cell classes are relevant for this stimulus direction,
    % as well as the parameters of the cortical, post-receptoral channel
    % second order filter
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

    nCellClasses = length(cellClasses);

    % Loop over eccentricities
    parfor ee=1:nEccs

        for cc = 1:nCellClasses

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

            % Grab the V1 MRI parameters that are fixed across
            % eccentricity
            secondOrderFc = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + 4);
            secondOrderQ = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + 5);

            % Grab the V1 MRI parameters that vary with eccentricity 
            v1SurroundDelay = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + nFixedParams + nEccs*0 + ee);
            v1SurroundIndex = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + nFixedParams + nEccs*1 + ee);
            v1Gain = pMRI((paramBlockIndex(cc)-1)*nParamsPerCellBlock + nFixedParams + nEccs*2 + ee);

            % Assemble the staged parameters
            surroundDelay = [lgnSurroundDelay v1SurroundDelay];
            surroundIndex = [lgnSurroundIndex v1SurroundIndex];
            gain = [lgnGain v1Gain];

            % Get the TTF
            thisResponse = ...
                returnV1TTFForEcc(cellClasses{cc},stimulusDirections{ss},...
                rgcTemporalModel,studiedEccentricites(ee),studiedFreqs,...
                stimulusContrastScale,surroundDelay,surroundIndex,gain,...
                secondOrderFc,secondOrderQ,nSubtractions);

            % Store this response
            responseMat(ss,ee,cc,:) = thisResponse;

        end
    end
end

response = [];
for ss = 1:nStims
    for ee=1:nEccs
        response = [response sum(squeeze(responseMat(ss,ee,:,:)))];
    end
end

end