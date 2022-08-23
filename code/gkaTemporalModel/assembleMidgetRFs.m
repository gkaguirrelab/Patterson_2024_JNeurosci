function [rfLMCone, rfBipolar, rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(cfCone, coneDelay, g, k, cfLowPass, cfInhibit, cf2ndStage, Q, midgetSurroundWeight, midgetSurroundDelay, chromaticCenterWeight, chromaticSurroundWeight)

syms f

rfLMCone = 10^g * stageFirstOrderLP(f,cfCone,1) .* stageFirstOrderLP(f,cfLowPass,1) .* stageDelay(f,coneDelay/1000);
rfBipolar = rfLMCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfMidgetLum = rfBipolar - midgetSurroundWeight.* rfBipolar .* stageDelay(f,midgetSurroundDelay/1000);
rfMidgetChrom = chromaticCenterWeight .* rfBipolar - ...
    midgetSurroundWeight .* chromaticSurroundWeight .* rfBipolar .* stageDelay(f,midgetSurroundDelay/1000);

end