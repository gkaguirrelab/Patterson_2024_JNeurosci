function [v1ChromAmplitude,v1ChromPhase] = returnV1ChromEccTTFFit(p,freqs,eccDeg)

% Unpack model parameters
secondOrderFc = p(1);
secondOrderQ = p(2);
surroundDelay = p(3);
surroundIndex = p(4);
gain = p(5);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1

% Obtain the response
rfChromV1Ecc = returnPostRetinalResponses(eccDeg,secondOrderFc,secondOrderQ,surroundIndex,surroundDelay,nSubtractions);
ttfComplex = double(subs(rfChromV1Ecc,freqs));
v1ChromAmplitude = abs(ttfComplex).*(1/length(rfChromV1Ecc));
v1ChromPhase = unwrap(angle(ttfComplex)).*(1/length(rfChromV1Ecc));

v1ChromAmplitude = gain.*v1ChromAmplitude;

end