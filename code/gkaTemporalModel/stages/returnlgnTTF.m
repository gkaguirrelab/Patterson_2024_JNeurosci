function [lgnAmplitude,lgnPhase] = returnlgnTTF(cellClass,stimulusDirection,rgcTemporalModel,surroundDelay,surroundIndex,gain,freqs)


% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:5:81;

% One surround subtraction at the LGN
nSubtractions = 1;

% Set up the call to the post retinal RF generator as a function of
% eccentricity
rfPostRetinalFunc = @(eccDeg) returnPostRetinalRF(...
    cellClass,stimulusDirection,rgcTemporalModel,eccDeg,nSubtractions,...
    surroundDelay,surroundIndex,gain,[],[]);

% Obtain the average TTF across the eccentricities
lgnAmplitude = zeros(size(freqs)); lgnPhase = zeros(size(freqs));
parfor ee=1:length(eccDegVals)

    % Get the lgn RF for this eccenticity
    rfPostRetinal = rfPostRetinalFunc(eccDegVals(ee));
    
    % Extract the amplitude and phase
    ttfComplex = double(subs(rfPostRetinal,freqs));
    lgnAmplitude = lgnAmplitude + abs(ttfComplex).*(1/length(eccDegVals));
    lgnPhase = lgnPhase + unwrap(angle(ttfComplex)).*(1/length(eccDegVals));
end


end