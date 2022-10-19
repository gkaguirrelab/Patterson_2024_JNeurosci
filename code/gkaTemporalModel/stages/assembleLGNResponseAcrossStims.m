function response = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*2;
LMRatio = pMRI(1);

% Loop through the stimulus directions and assemble the response
responseMat = zeros(nStims,nFreqs);

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

    responseAtStim = zeros(1,length(studiedFreqs));
    for cc = 1:length(cellClasses)

        % What is the index in the param structure of this cell class?
        cellClassIndex = find(strcmp(cellClassOrder,cellClasses{cc}));

        % Grab the block of pMRI parameters that correspond to this cell class
        pMRICellBlock = pMRI( 1+nUniqueParams+(cellClassIndex-1)*nParamsPerCellBlock:nUniqueParams+cellClassIndex*nParamsPerCellBlock);

        % Get the TTF
        responseAtStim = responseAtStim + ...
            returnlgnTTF(cellClasses{cc},stimulusDirections{ss},...
            rgcTemporalModel,pMRICellBlock,LMRatio,...
            studiedFreqs,studiedEccentricites,nFixedParams);

    end

    % Store this response
    responseMat(ss,:) = responseAtStim ./ length(cellClasses);

end

response = [];
for ss = 1:nStims
    response = [response squeeze(responseMat(ss,:))];
end

end
