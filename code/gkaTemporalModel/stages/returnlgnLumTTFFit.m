function [lgnLumAmplitude,lgnLumPhase] = returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 1; % One delayed surround stage at the LGN

% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:20:61;

% Interpolate the surround index values for the modeled eccentricities
surroundIndex = interp1(log10(unique(v1Eccentricity)),surroundIndex,log10(eccDegVals),'linear','extrap');

% Loop through eccentricities and obtain modeled responses
rfLumLGN = {};
parfor ee=1:length(eccDegVals)
    [~,rfLumLGN{ee}] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
end
rfLumLGN = reshape(vertcat(rfLumLGN{:}),length(eccDegVals)*2,1);

% Obtain the average luminance TTF across these eccentricities
lgnLumAmplitude = zeros(size(lgnFreqX)); lgnLumPhase = zeros(size(lgnFreqX));
for ii=1:length(rfLumLGN)
    ttfComplex = double(subs(rfLumLGN{ii},lgnFreqX));
    lgnLumAmplitude = lgnLumAmplitude + abs(ttfComplex).*(1/length(rfLumLGN));
    lgnLumPhase = lgnLumPhase + unwrap(angle(ttfComplex)).*(1/length(rfLumLGN));
end

lgnLumAmplitude = lgnLumAmplitude.*lgnGain;

end