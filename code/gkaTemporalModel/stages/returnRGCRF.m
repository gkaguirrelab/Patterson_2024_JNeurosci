function [rfRGC, rfBipolar, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight)
% returnRGCRF
%
% Derive and return the complex Fourier domain symbolic equations that
% express temporal sensitiviy for a retinal ganglion cell

% Extract the parameters
g = pBlock(1); k = pBlock(2);
cfInhibit = pBlock(3); cf2ndStage = pBlock(4); Q = pBlock(5);
surroundWeight = pBlock(6); surroundDelay = pBlock(7);

% Assemble the Fourier domain representation
syms f
rfCone = 10^4 * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfBipolar = rfCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfRGC = chromaticCenterWeight .* rfBipolar - ...
    surroundWeight .* chromaticSurroundWeight .* rfBipolar .* stageDelay(f,surroundDelay/1000);
rfRGC = g.*rfRGC;

end
