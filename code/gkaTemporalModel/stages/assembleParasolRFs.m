function [rfLMCone, rfParasolBipolar, rfParasolLum, f] = assembleParasolRFs(cfCone, coneDelay, g, k, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay, f)

if nargin<10
    syms f
end

rfLMCone = 10^g * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfParasolBipolar = rfLMCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfParasolLum = rfParasolBipolar - surroundWeight.* rfParasolBipolar .* stageDelay(f,surroundDelay/1000);

end