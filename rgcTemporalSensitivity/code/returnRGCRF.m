function [rfRGC, rfBipolar, rfCone] = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight)
% returnRGCRF
%
% Derive and return the complex Fourier domain symbolic equations that
% express temporal sensitiviy for a retinal ganglion cell
%
% Examples:
%{
    pRGCBlock = [ 0.2469    0.5922    5.8057   39.0848    1.3149    0.8984    2.1112 ];
    cfCone = 14.2376;
    coneDelay = 13.7579;
    chromaticCenterWeight = 0.9960;
    chromaticSurroundWeight = 0.1410;
    rfRGC = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);
    plotRF(rfRGC);
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
