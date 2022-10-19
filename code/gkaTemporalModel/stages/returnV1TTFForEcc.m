function [v1Amplitude,v1Phase] = returnV1TTFForEcc(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,pMRIBlock,LMRatio,freqs)

% Unpack model parameters
v1Gain = pMRIBlock(end);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1

% Get the post retinal temporal RF model
rfPostRetinal = returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,pMRIBlock,LMRatio,nSubtractions);

% Derive the amplitude and phase from the Fourier model, combining the two
% LMS components if present
ttfComplex = double(subs(rfPostRetinal,freqs));
v1Amplitude = abs(ttfComplex);
v1Phase = unwrap(angle(ttfComplex));

% Apply the gain
v1Amplitude = v1Gain.*v1Amplitude;

end