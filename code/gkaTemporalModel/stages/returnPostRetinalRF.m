function [rfPostRetinal, rfRGC, rfBipolar, rfCone] = returnPostRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,pMRI,nSubtractions)

% Unpack the MRI search parameters
LMConeRatio = pMRI(1);
% lgnGain = pMRI(2); unused
secondOrderFc = pMRI(3);
secondOrderQ = pMRI(4);
surroundDelay = pMRI(5);
surroundIndex = pMRI(6);
% v1Gain = pMRI(7); used by the calling function

% Obtain the chromatic weights. Note that the LMConeRatio that is passed
% by the pMRI model parameters overrides the LMConeRatio that was used to
% define the pRGC model.
switch cellClass
    case 'midget'
        switch stimulusDirection
            case 'LminusM'
                [chromaticCenterWeight,chromaticSurroundWeight] = ...
                    returnMidgetChromaticWeights(eccDeg,LMConeRatio);
            case 'LMS'
                chromaticCenterWeight = 1; chromaticSurroundWeight = 1;
            case 'S'
                chromaticCenterWeight = 0; chromaticSurroundWeight = 0;
        end
    case 'parasol'
        switch stimulusDirection
            case 'LminusM'
                chromaticCenterWeight = 0; chromaticSurroundWeight = 0;
            case 'LMS'
                chromaticCenterWeight = 1; chromaticSurroundWeight = 1;
            case 'S'
                chromaticCenterWeight = 0; chromaticSurroundWeight = 0;
        end
    case 'bistratified'
        switch stimulusDirection
            case 'LminusM'
                chromaticCenterWeight = 0; chromaticSurroundWeight = 0;
            case 'LMS'
                chromaticCenterWeight = 0; chromaticSurroundWeight = 0;
            case 'S'
                chromaticCenterWeight = 1; chromaticSurroundWeight = 1;
        end
end

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
