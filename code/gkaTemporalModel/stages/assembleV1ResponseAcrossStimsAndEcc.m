function [response, responseMat] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nUniqueParams,nFixedParams)

% Identify the studied eccentricities and stimulus frequencies
nStims = length(stimulusDirections);
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nParamsPerCellBlock = nFixedParams+nEccs*2;
nSubtractions = 2; % Perform two levels of surround-delayed subtraction for V1

% The shared LGN parameters
lgnSurroundDelay = pMRI(1);
lgnSurroundIndex = pMRI(2);

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
            cellClassIndex = 1;
        case 'S'
            cellClasses = {'bistratified'};
            cellClassIndex = 2;
        case 'LMS'
            cellClasses = {'parasol','midget'};
            cellClassIndex = [3 4];
    end

    nCellClasses = length(cellClasses);

    % Loop over eccentricities
    parfor ee=1:nEccs

        for cc = 1:nCellClasses

            % Grab the LGN parameters, which are organized by RGC class
            lgnGain = [];
            switch cellClasses{cc}
                case 'midget'
                    lgnGain = pMRI(3);
                case 'bistratified'
                    lgnGain = pMRI(4);
                case 'parasol'
                    lgnGain = pMRI(5);
            end

            % Grab "fixed" MRI parameters for this post-receptoral path,
            % which vary by cell type and stimulus direction
            secondOrderFc = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+1);
            secondOrderQ = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+2);
            v1SurroundDelay = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+3);

            % Grab the "floating" parameters that vary with eccentricity
            v1SurroundIndex = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+ee);
            v1Gain = pMRI(nUniqueParams+(cellClassIndex(cc)-1)*nParamsPerCellBlock+nFixedParams+nEccs+ee);

            % Assemble the staged parameters
            surroundDelay = [lgnSurroundDelay v1SurroundDelay];
            surroundIndex = [lgnSurroundIndex v1SurroundIndex];
            gain = [lgnGain v1Gain];

            % Get the TTF
            thisResponse = ...
                returnV1TTFForEcc(cellClasses{cc},stimulusDirections{ss},...
                rgcTemporalModel,studiedEccentricites(ee),studiedFreqs,...
                surroundDelay,surroundIndex,gain,secondOrderFc,secondOrderQ,nSubtractions);

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