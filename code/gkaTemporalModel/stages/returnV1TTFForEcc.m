function [v1Amplitude,v1Phase] = returnV1TTFForEcc(stimulusDirection,rgcTemporalModel,eccDeg,pMRI,freqs)

% Unpack model parameters
v1Gain = pMRI(7);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1

% Get the post retinal temporal RF models
rfPostRetinal = {};
switch stimulusDirection
    case 'LminusM'
        rfPostRetinal{1} = returnPostRetinalRF('midget',stimulusDirection,rgcTemporalModel,eccDeg,pMRI,nSubtractions);
    case 'S'
        rfPostRetinal{1} = returnPostRetinalRF('bistratified',stimulusDirection,rgcTemporalModel,eccDeg,pMRI,nSubtractions);
    case 'LMS'
        rfPostRetinal{1} = returnPostRetinalRF('midget',stimulusDirection,rgcTemporalModel,eccDeg,pMRI,nSubtractions);
        rfPostRetinal{2} = returnPostRetinalRF('parasol',stimulusDirection,rgcTemporalModel,eccDeg,pMRI,nSubtractions);
end

% Derive the amplitude and phase from the Fourier model, combining the two
% LMS components if present
v1Amplitude=zeros(size(freqs)); v1Phase=zeros(size(freqs));
for ii=1:length(rfPostRetinal)
    ttfComplex = double(subs(rfPostRetinal{ii},freqs));
    v1Amplitude = v1Amplitude + abs(ttfComplex).*(1/length(rfPostRetinal));
    v1Phase = v1Phase + unwrap(angle(ttfComplex)).*(1/length(rfPostRetinal));
end

% Apply the gain
v1Amplitude = v1Gain.*v1Amplitude;

end