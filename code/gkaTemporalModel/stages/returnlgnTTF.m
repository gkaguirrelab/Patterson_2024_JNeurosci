function [lgnAmplitude,lgnPhase] = returnlgnTTF(stimulusDirection,rgcTemporalModel,pMRI,lgnFreqX,v1Eccentricity)

% Some things we need to know
nFixed = 5;
nEcc = length(unique(v1Eccentricity));

% Unpack model parameters
lgnGain = pMRI(2);
surroundIndexV1 = pMRI(nFixed+1:nFixed+nEcc);
nSubtractions = 1; % One delayed surround stage at the LGN

% An anonymous function that returns the surround index as a function of
% eccentriciy, interpolated from the values measured in V1
interpSurroundIndex = @(eccDeg) interp1(log10(unique(v1Eccentricity)),surroundIndexV1,log10(eccDeg),'linear','extrap');

% Anonymous function to assemble a pMRI vector in which the surround index
% has been replaced with the interpolated value
mypMRI = @(eccDeg) [pMRI(1:nFixed) interpSurroundIndex(eccDeg)];

% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:20:61;

% Set up the call to the post retinal RF generator as a function of
% eccentricity
rfPostRetinalFunc = {};
switch stimulusDirection
    case 'LminusM'
        rfPostRetinalFunc{1} = @(eccDeg) returnPostRetinalRF('midget',stimulusDirection,rgcTemporalModel,eccDeg,mypMRI(eccDeg),nSubtractions);
    case 'S'
        rfPostRetinalFunc{1} = @(eccDeg) returnPostRetinalRF('bistratified',stimulusDirection,rgcTemporalModel,eccDeg,mypMRI(eccDeg),nSubtractions);
    case 'LMS'
        rfPostRetinalFunc{1} = @(eccDeg) returnPostRetinalRF('midget',stimulusDirection,rgcTemporalModel,eccDeg,mypMRI(eccDeg),nSubtractions);
        rfPostRetinalFunc{2} = @(eccDeg) returnPostRetinalRF('parasol',stimulusDirection,rgcTemporalModel,eccDeg,mypMRI(eccDeg),nSubtractions);
end

% Obtain the average luminance TTF across these eccentricities
lgnAmplitude = zeros(size(lgnFreqX)); lgnPhase = zeros(size(lgnFreqX));
for ii=1:length(rfPostRetinalFunc)
    parfor ee=1:length(eccDegVals)
        % Get the lgn RF for this eccenticity
        rfPostRetinal = rfPostRetinalFunc{ii}(eccDegVals(ee));
        % Extract the amplitude and phase
        ttfComplex = double(subs(rfPostRetinal,lgnFreqX));
        lgnAmplitude = lgnAmplitude + abs(ttfComplex).*(1/length(rfPostRetinalFunc));
        lgnPhase = lgnPhase + unwrap(angle(ttfComplex)).*(1/length(rfPostRetinalFunc));
    end
end

% Apply the LGN gain parameters
lgnAmplitude = lgnAmplitude.*lgnGain;

end