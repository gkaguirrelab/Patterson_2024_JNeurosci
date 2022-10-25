function [response, responseMat] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nUniqueParams,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*2;
nSubtractions = 2; % Perform two levels of surround-delayed subtraction for V1

responseMat = zeros(nStims,nEccs,2,nFreqs);

% Loop through the stimulus directions and assemble the response
for ss = 1:nStims

    % Identify which cell classes are relevant for this stimulus direction,
    % as well as the parameters of the cortical, post-receptoral channel
    % second order filter
    switch stimulusDirections{ss}
        case 'LminusM'
            cellClasses = {'midget'};
            cellClassIndex = 1;
            stimulusIndex = 1;
        case 'S'
            cellClasses = {'bistratified'};
            cellClassIndex = 2;
            stimulusIndex = 2;
        case 'LMS'
            cellClasses = {'parasol','midget'};
            cellClassIndex = [3 4];
            stimulusIndex = 3;
    end

    secondOrderFc = pMRI( (stimulusIndex-1)*(nUniqueParams/3)+1 );
    secondOrderQ = pMRI( (stimulusIndex-1)*(nUniqueParams/3)+2 );
    surroundDelay = pMRI( (stimulusIndex-1)*(nUniqueParams/3)+3 );

    nCellClasses = length(cellClasses);

    % Loop over eccentricities
    parfor ee=1:nEccs

        for cc = 1:nCellClasses

            % Grab MRI parameters that correspond to this cell class,
            % stimulus, and eccentricity
            surroundIndex = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+ee);
            v1Gain = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+nEccs+ee);

            % Get the TTF
            thisResponse = ...
                returnV1TTFForEcc(cellClasses{cc},stimulusDirections{ss},...
                rgcTemporalModel,studiedEccentricites(ee),studiedFreqs,...
                v1Gain,surroundDelay,surroundIndex,secondOrderFc,secondOrderQ,nSubtractions);

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