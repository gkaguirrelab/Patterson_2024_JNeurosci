function [rfPostRetinal, rfRGC, rfBipolar, rfCone] = returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,pMRI,LMRatio,nSubtractions)

% Unpack the MRI search parameters
secondOrderFc = pMRI(2);
secondOrderQ = pMRI(3);
surroundDelay = pMRI(4);
surroundIndex = pMRI(5);

% Obtain the chromatic weights
[chromaticCenterWeight,chromaticSurroundWeight] = ...
    returnRGCChromaticWeights(cellClass,stimulusDirection,eccDeg,LMRatio);

% Obtain the pRGC model params for this cell class and eccentricity
cellIndex = strcmp(cellClass,{'midget','parasol','bistratified'});
pBlockRGC = nan(1,7);
for ii = 1:7
    pBlockRGC(ii) = rgcTemporalModel.pFitByEccen{ii,cellIndex}(eccDeg);
end

% Get the temporal receptive field for this RGC
[rfRGC, rfBipolar, rfCone] = ...
    returnRGCRF(pBlockRGC,rgcTemporalModel.cfCone,rgcTemporalModel.coneDelay,chromaticCenterWeight,chromaticSurroundWeight);

% Copy the RGC model into the post-retinal variables
rfPostRetinal = rfRGC;

% Apply iterative, delayed surround subtraction. The expectation is one
% iteration for LGN, two iterations for V1
syms f
for ii = 1:nSubtractions
    rfPostRetinal = rfPostRetinal - surroundIndex * rfPostRetinal.*stageDelay(f,surroundDelay/1000);
end

% Apply a post-lgn, second order low pass filter
if nSubtractions > 1
    rfPostRetinal = rfPostRetinal.*stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
end

% Drasdo 2007 equation for the midget fraction as a function of
% eccentricity
midgetFraction = @(eccDeg) 0.8928*(1+eccDeg/41.03).^(-1);

% At each eccentricity, what fraction of the total number of RGCs is of a
% given cell class? We simply assume that the bistratified RGCs make up a
% uniform 10% of all RGCs across eccentricity
switch cellClass
    case 'midget'
        proportionFunc = @(ecc) midgetFraction(ecc);
    case 'parasol'
        proportionFunc = @(ecc) 0.9-midgetFraction(ecc);
    case 'bistratified'
        proportionFunc = @(ecc) 0.1;
end

% Loop over the meridians and obtain the RGCf density functions
meridianAngles = [0 90 180 270];
for mm = 1:4
    totalRGCfDensity{mm} = watsonTotalRFDensityByEccenDegVisual(meridianAngles(mm));
end

% Obtain the mean, total ganglion receptive field density as a function of
% eccentricity
totalRGCfDensityAtEcc = @(eccDeg) mean(cellfun(@(x) x(eccDeg),totalRGCfDensity));

% Scale by cell density at this eccentricity, accounting for the
% change in surface area as a function of eccentricity (eccDeg/90), and
% the fraction of cells of each class.
rfPostRetinal = rfPostRetinal ...
    .* totalRGCfDensityAtEcc(eccDeg) .* (eccDeg/90) ...
    .* proportionFunc(eccDeg);

end
