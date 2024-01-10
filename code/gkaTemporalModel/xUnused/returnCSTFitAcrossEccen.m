function [response,rfsAtEcc] = returnCSTFitAcrossEccen(p,stimulusDirections,studiedEccentricites,studiedFreqs)
% Assemble the response across eccentricity locations

% Fixed params of the analysis
nCells = 3;
nParams = 3;

% Loop over the passed eccentricities
parfor ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    idx = [1,(ee-1)*nParams*nCells+2 : ee*nParams*nCells+1];
    subP = p(idx);

    % Obtain the response at this eccentricity
    [thisTTF,rfsAtEcc{ee}] = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs);

    % Detect if the response is not band pass and in that case make it a
    % bad fit so that we avoid finding these solutions
    for ss = 1:length(stimulusDirections)
        [~,idx] = max(thisTTF(ss,:));
        if idx < length(studiedFreqs)/3
            thisTTF(ss,:) = 100;
        end
    end

    % Store this loop result
    ttfAtEcc{ee} = thisTTF;

end

% Pull the fits out of the par pool cell array and place in a matrix
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = ttfAtEcc{ee};
end

end