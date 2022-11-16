function [rfRGC, rfBipolar, rfCone] = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight)
% returnRGCRF
%
% Derive and return the complex Fourier domain symbolic equations that
% express temporal sensitiviy for a retinal ganglion cell
%
% Examples:
%{
    [rfRGC, rfBipolar, rfCone] = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight)
%}

% Extract the parameters
g = pRGCBlock(1); k = pRGCBlock(2);
cfInhibit = pRGCBlock(3); cf2ndStage = pRGCBlock(4); Q = pRGCBlock(5);
surroundWeight = pRGCBlock(6); surroundDelay = pRGCBlock(7);

% Assemble the Fourier domain representation
rfCone = 10^4 * stageFirstOrderLP(cfCone,2) .* stageDelay(coneDelay/1000);
rfBipolar = rfCone .* stageInhibit(cfInhibit,k) .* stageSecondOrderLP(cf2ndStage,Q);
rfRGC = chromaticCenterWeight .* rfBipolar - ...
    surroundWeight .* chromaticSurroundWeight .* rfBipolar .* stageDelay(surroundDelay/1000);
rfRGC = g.*rfRGC;

end
