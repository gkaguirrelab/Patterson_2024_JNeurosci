function response = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nUniqueParams,nFixedParams)

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
            cellClassIndex = 1;
        case 'S'
            cellClasses = {'bistratified'};
            cellClassIndex = 2;
        case 'LMS'
            cellClasses = {'parasol','midget'};
            cellClassIndex = [3 4];
    end

    % Loop over eccentricities
    parfor ee=1:nEccs

        % A variable to hold the response
        responseAtEcc = zeros(1,nFreqs);

        for cc = 1:length(cellClasses)

            % Grab MRI parameters that correspond to this cell class,
            % stimulus, and eccentricity
            secondOrderFc = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+2)
            secondOrderQ = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+3)
            surroundDelay = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+4)
            surroundIndex = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+4+ee);
            v1Gain = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+4+nEccs+ee);

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