function [rfRGC, rfBipolar, rfCone] = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight)
% returnRGCRF
%
% Derive and return the complex Fourier domain symbolic equations that
% express temporal sensitiviy for a retinal ganglion cell

% Extract the parameters
g = pRGCBlock(1); k = pRGCBlock(2);
cfInhibit = pRGCBlock(3); cf2ndStage = pRGCBlock(4); Q = pRGCBlock(5);
surroundWeight = pRGCBlock(6); surroundDelay = pRGCBlock(7);

% Assemble the Fourier domain representation
syms f
rfCone = 10^4 * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfBipolar = rfCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfRGC = chromaticCenterWeight .* rfBipolar - ...
    surroundWeight .* chromaticSurroundWeight .* rfBipolar .* stageDelay(f,surroundDelay/1000);
rfRGC = g.*rfRGC;

end
