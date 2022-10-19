function response = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*2;
LMRatio = pMRI(1);

% Loop through the stimulus directions and assemble the response
responseMat = zeros(nStims,nEccs,nFreqs);

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

    % Loop over eccentricities
    parfor ee=1:nEccs

        % A variable to hold the response
        responseAtEcc = zeros(1,nFreqs);

        for cc = 1:length(cellClasses)

            % What is the index in the param structure of this cell class?
            cellClassIndex = find(strcmp(cellClassOrder,cellClasses{cc}));

            % Grab the block of pMRI parameters that correspond to this
            % cell class
            pMRICellBlock = pMRI( 1+nUniqueParams+(cellClassIndex-1)*nParamsPerCellBlock:nUniqueParams+cellClassIndex*nParamsPerCellBlock);

            % Grab the portion of pMRI cell that corresponds to this
            % eccentriciy
            pMRICellEccBlock = [pMRICellBlock(1:nFixedParams) pMRICellBlock(nFixedParams+ee) pMRICellBlock(nFixedParams+nEccs+ee)];

            % Get the TTF
            responseAtEcc = responseAtEcc + ...
                returnV1TTFForEcc(cellClasses{cc},stimulusDirections{ss},...
                rgcTemporalModel,studiedEccentricites(ee),pMRICellEccBlock,LMRatio,studiedFreqs);
        end

        % Store this response
        responseMat(ss,ee,:) = responseAtEcc ./ length(cellClasses);

    end
end

response = [];
for ss = 1:nStims
    for ee=1:nEccs
        response = [response squeeze(responseMat(ss,ee,:))'];
    end
end

end