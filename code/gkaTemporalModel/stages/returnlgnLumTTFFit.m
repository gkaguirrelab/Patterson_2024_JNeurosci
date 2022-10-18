function [lgnLumAmplitude,lgnLumPhase] = returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity)

nFixed = 5;
nEcc = 6;

% Unpack model parameters
LMConeRatio = p(1);
lgnGain = p(2);
secondOrderFc = p(3);
secondOrderQ = p(4);
surroundDelay = p(5);
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
    [~,rfLumLGN{ee}] = returnPostRetinalResponses(eccDegVals(ee),LMConeRatio,secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
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