function [v1Amplitude,v1Phase] = returnV1TTFForEcc(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,freqs,stimulusContrastScale,surroundDelay,surroundIndex,gain,secondOrderFc,secondOrderQ,nSubtractions)

% Get the post retinal temporal RF model
rfPostRetinal = returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,nSubtractions,stimulusContrastScale,surroundDelay,surroundIndex,gain,secondOrderFc,secondOrderQ);

% Derive the amplitude and phase from the Fourier model
ttfComplex = double(subs(rfPostRetinal,freqs));
v1Amplitude = abs(ttfComplex);
v1Phase = unwrap(angle(ttfComplex));

end