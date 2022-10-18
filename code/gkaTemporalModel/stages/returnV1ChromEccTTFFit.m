function [v1ChromAmplitude,v1ChromPhase] = returnV1ChromEccTTFFit(p,freqs,eccDeg)

% Unpack model parameters
LMConeRatio = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(5);
gain = p(6);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1

% Obtain the response
rfChromV1Ecc = returnPostRetinalResponses(eccDeg,LMConeRatio,secondOrderFc,secondOrderQ,surroundIndex,surroundDelay,nSubtractions);
ttfComplex = double(subs(rfChromV1Ecc,freqs));
v1ChromAmplitude = abs(ttfComplex).*(1/length(rfChromV1Ecc));
v1ChromPhase = unwrap(angle(ttfComplex)).*(1/length(rfChromV1Ecc));

v1ChromAmplitude = gain.*v1ChromAmplitude;

end