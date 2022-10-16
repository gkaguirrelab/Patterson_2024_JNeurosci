function [rfMidgetChrom, rfMidgetLum, rfLMCone, rfMidgetBipolar, f] = parseParamsMidget(pBlock, cfCone, coneDelay, LMRatio, eccDeg)
% parseParamsMidget
% Given the parameters for an eccentricity location, as well as the fixed
% cone-stage parameters and the LM ratio, derive and return the complex
% Fourier domain symbolic equations the express temporal sensitiviy of the
% midgets.

% Extract the parameters
g = pBlock(1); k = pBlock(2);
cfInhibit = pBlock(3); cf2ndStage = pBlock(4); Q = pBlock(5);
surroundWeight = pBlock(6); surroundDelay = pBlock(7);

% Simulate 1000 RGCs with random sampling from the cone mosaic at the
% specified eccentricity and with the specified LMRatio. For each RGC,
% obtain the chromatic weight on the center and surround. We then take the
% absolute value of these, as we assume that the behavior of an L-center
% RGC is the same as the behavior of an M-center RGC.
tmpCenterWeight = [];
tmpSurroundWeight = [];
for cc = 1:1000
    [tmpC,tmpS] = ...
        returnChromaticWeights(eccDeg,LMRatio);
    tmpCenterWeight(cc) = tmpC;
    tmpSurroundWeight(cc) = tmpS;
end
chromaticCenterWeight = mean(abs(tmpCenterWeight));
chromaticSurroundWeight = mean(abs(tmpSurroundWeight));

% Assemble the temporal receptive fields
[rfLMCone, rfMidgetBipolar, rfMidgetChrom, rfMidgetLum, f] = assembleMidgetRFs(...
    cfCone, coneDelay, ...
    g, k, cfInhibit, cf2ndStage, Q, ...
    surroundWeight, surroundDelay, ...
    chromaticCenterWeight, chromaticSurroundWeight);

end
