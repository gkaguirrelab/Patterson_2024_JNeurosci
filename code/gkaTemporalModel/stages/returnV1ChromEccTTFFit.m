function [v1ChromAmplitude,v1ChromPhase] = returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1
v1LumGain = p(nFixed+nEcc+1:end);
v1ChromGain = v1LumGain;

% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
v1ChromAmplitude=[];v1ChromPhase=[];

parfor ee=1:length(eccDegVals)
    [rfChromV1Ecc] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
    ttfComplex = double(subs(rfChromV1Ecc,studiedFreqs));
    v1ChromAmplitude(ee,:) = v1ChromGain(ee).*abs(ttfComplex);
    v1ChromPhase(ee,:) = unwrap(angle(ttfComplex));
end

v1ChromAmplitude = reshape(v1ChromAmplitude',1,nEcc*length(studiedFreqs));
v1ChromPhase = reshape(v1ChromPhase',1,nEcc*length(studiedFreqs));

end