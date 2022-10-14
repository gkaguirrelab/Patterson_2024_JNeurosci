function [rfLMCone, rfBipolar, rfParasaolLum, f] = assembleParasolRFs(cfCone, coneDelay, g, k, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay, f)

if nargin<10
    syms f
end

rfLMCone = 10^g * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfBipolar = rfLMCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfParasaolLum = rfBipolar - surroundWeight.* rfBipolar .* stageDelay(f,surroundDelay/1000);

end