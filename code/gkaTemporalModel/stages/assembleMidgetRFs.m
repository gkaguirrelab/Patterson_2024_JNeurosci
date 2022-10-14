function [rfLMCone, rfBipolar, rfMidgetChrom, rfMidgetLum, f] = assembleMidgetRFs(cfCone, coneDelay, g, k, cfInhibit, cf2ndStage, Q, midgetSurroundWeight, midgetSurroundDelay, chromaticCenterWeight, chromaticSurroundWeight,f)

if nargin<12
    syms f
end

rfLMCone = 10^g * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfBipolar = rfLMCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfMidgetLum = rfBipolar - midgetSurroundWeight.* rfBipolar .* stageDelay(f,midgetSurroundDelay/1000);
rfMidgetChrom = chromaticCenterWeight .* rfBipolar - ...
    midgetSurroundWeight .* chromaticSurroundWeight .* rfBipolar .* stageDelay(f,midgetSurroundDelay/1000);

end