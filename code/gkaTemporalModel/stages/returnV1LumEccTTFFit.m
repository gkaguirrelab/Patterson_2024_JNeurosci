function [v1LumAmplitude,v1LumPhase] = returnV1LumEccTTFFit(p,freqs,eccDeg)

% Unpack model parameters
LMConeRatio = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(5);
gain = p(6);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1

v1LumAmplitude=zeros(size(freqs)); v1LumPhase=zeros(size(freqs));

[~,rfLumV1Ecc] = returnPostRetinalResponses(eccDeg,LMConeRatio,secondOrderFc,secondOrderQ,surroundIndex,surroundDelay,nSubtractions);
for ii=1:length(rfLumV1Ecc)
    ttfComplex = double(subs(rfLumV1Ecc{ii},freqs));
    v1LumAmplitude = v1LumAmplitude + abs(ttfComplex).*(1/length(rfLumV1Ecc));
    v1LumPhase = v1LumPhase + unwrap(angle(ttfComplex)).*(1/length(rfLumV1Ecc));
end

v1LumAmplitude = gain.*v1LumAmplitude;

end