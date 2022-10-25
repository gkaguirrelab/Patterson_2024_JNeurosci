function [lgnAmplitude,lgnPhase] = returnlgnTTF(cellClass,stimulusDirection,rgcTemporalModel,pMRICellBlock,studiedFreqs,studiedEccentricites,nFixed)

% Some things we need to know
nV1Eccs = length(studiedEccentricites);

% Unpack model parameters
lgnGain = pMRICellBlock(1);
surroundDelay = pMRICellBlock(4);
surroundIndexV1 = pMRICellBlock(nFixed+1:nFixed+nV1Eccs);
nSubtractions = 1; % One delayed surround stage at the LGN

% An anonymous function that returns the surround index as a function of
% eccentriciy, interpolated from the values measured in V1
interpSurroundIndex = @(eccDeg) interp1(log10(studiedEccentricites),surroundIndexV1,log10(eccDeg),'linear','extrap');

% Anonymous function to assemble a pMRI vector in which the surround index
% has been replaced with the interpolated value
mypMRI = @(eccDeg) [pMRICellBlock(1:nFixed) interpSurroundIndex(eccDeg)];

% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:5:81;

% Set up the call to the post retinal RF generator as a function of
% eccentricity
rfPostRetinalFunc = @(eccDeg) returnPostRetinalRF(...
    cellClass,stimulusDirection,rgcTemporalModel,eccDeg,nSubtractions,...
    surroundDelay,interpSurroundIndex(eccDeg),[],[]);

% Obtain the average TTF across the eccentricities
lgnAmplitude = zeros(size(studiedFreqs)); lgnPhase = zeros(size(studiedFreqs));
parfor ee=1:length(eccDegVals)
    % Get the lgn RF for this eccenticity
    rfPostRetinal = rfPostRetinalFunc(eccDegVals(ee));
    % Extract the amplitude and phase
    ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
    lgnAmplitude = lgnAmplitude + abs(ttfComplex).*(1/length(eccDegVals));
    lgnPhase = lgnPhase + unwrap(angle(ttfComplex)).*(1/length(eccDegVals));
end

% Apply the LGN gain parameters
lgnAmplitude = lgnAmplitude.*lgnGain;

end