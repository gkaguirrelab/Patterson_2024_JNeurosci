%% parseParamsParasol
function [rfParasaolLum, rfLMCone, rfBipolar, f] = parseParamsParasol(pBlock, cfCone, coneDelay)
% Given the parameters for an eccentricity location, as well as the fixed
% cone-stage parameters, derive and return the complex
% Fourier domain symbolic equations the express temporal sensitiviy of the
% midgets.

% Extract the parameters
g = pBlock(1); k = pBlock(2);
cfInhibit = pBlock(3); cf2ndStage = pBlock(4); Q = pBlock(5);
surroundWeight = pBlock(6); surroundDelay = pBlock(7);

% Assemble the temporal receptive fields
[rfLMCone, rfBipolar, rfParasaolLum, f] = assembleParasolRFs(...
    cfCone, coneDelay, ...
    g, k, cfInhibit, cf2ndStage, Q, ...
    surroundWeight, surroundDelay);

end