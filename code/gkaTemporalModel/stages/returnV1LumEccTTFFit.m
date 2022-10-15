function [v1LumAmplitude,v1LumPhase] = returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity)

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

% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
v1LumAmplitude=[];v1LumPhase=[];

parfor ee=1:length(eccDegVals)
    [~,rfLumV1Ecc] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
    gainVals = zeros(size(studiedFreqs)); angleVals = zeros(size(studiedFreqs));
    for ii=1:length(rfLumV1Ecc)
        ttfComplex = double(subs(rfLumV1Ecc{ii},studiedFreqs));
        gainVals = gainVals + abs(ttfComplex).*(1/length(rfLumV1Ecc));
        angleVals = angleVals + unwrap(angle(ttfComplex)).*(1/length(rfLumV1Ecc));
    end
    v1LumAmplitude(ee,:) = v1LumGain(ee).*gainVals;
    v1LumPhase(ee,:) = angleVals;
end

v1LumAmplitude = reshape(v1LumAmplitude',1,nEcc*length(studiedFreqs));
v1LumPhase = reshape(v1LumPhase',1,nEcc*length(studiedFreqs));

end