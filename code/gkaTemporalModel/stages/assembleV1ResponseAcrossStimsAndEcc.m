function response = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*2;

% Loop through the stimulus directions and assemble the response
responseMat = zeros(nStims,nEccs,nFreqs);

for ss = 1:nStims

    % Identify which cell classes are relevant for this stimulus direction,
    % as well as the parameters of the cortical, post-receptoral channel
    % second order filter
    switch stimulusDirections{ss}
        case 'LminusM'
            cellClasses = {'midget'};
            secondOrderFc = pMRI(1); secondOrderQ = pMRI(2);
        case 'S'
            cellClasses = {'bistratified'};
            secondOrderFc = pMRI(3); secondOrderQ = pMRI(4);
        case 'LMS'
            cellClasses = {'parasol','midget'};
            secondOrderFc = pMRI(5); secondOrderQ = pMRI(6);
    end

    % Loop over eccentricities
    parfor ee=1:nEccs

        % A variable to hold the response
        responseAtEcc = zeros(1,nFreqs);

        for cc = 1:length(cellClasses)

            % What is the index in the param structure of this cell class?
            cellClassIndex = find(strcmp(cellClassOrder,cellClasses{cc}));

            % Grab the block of pMRI parameters that correspond to this
            % cell class and eccentricity
            surroundDelay = pMRI(1+nUniqueParams+(cellClassIndex-1)*nParamsPerCellBlock+1);
            surroundIndex = pMRI(1+nUniqueParams+(cellClassIndex-1)*nParamsPerCellBlock+1+ee);
            v1Gain = pMRI(1+nUniqueParams+(cellClassIndex-1)*nParamsPerCellBlock+nEccs+1+ee);

            % Get the TTF
            thisResponse = ...
                returnV1TTFForEcc(cellClasses{cc},stimulusDirections{ss},...
                rgcTemporalModel,studiedEccentricites(ee),studiedFreqs,...
                v1Gain,surroundDelay,surroundIndex,secondOrderFc,secondOrderQ);

            % Add this response to the total response
            responseAtEcc = responseAtEcc + thisResponse;
        end

        % Store this response
        responseMat(ss,ee,:) = responseAtEcc;

    end
end

response = [];
for ss = 1:nStims
    for ee=1:nEccs
        response = [response squeeze(responseMat(ss,ee,:))'];
    end
end

end