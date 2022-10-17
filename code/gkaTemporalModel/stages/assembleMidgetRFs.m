function [rfLMCone, rfMidgetBipolar, rfMidgetChrom, rfMidgetLum, f] = assembleMidgetRFs(cfCone, coneDelay, g, k, cfInhibit, cf2ndStage, Q, midgetSurroundWeight, midgetSurroundDelay, chromaticCenterWeight, chromaticSurroundWeight,f)

if nargin<12
    syms f
end

rfLMCone = 10^4 * stageFirstOrderLP(f,cfCone,2) .* stageDelay(f,coneDelay/1000);
rfMidgetBipolar = rfLMCone .* stageInhibit(f,cfInhibit,k) .* stageSecondOrderLP(f,cf2ndStage,Q);
rfMidgetLum = rfMidgetBipolar - midgetSurroundWeight.* rfMidgetBipolar .* stageDelay(f,midgetSurroundDelay/1000);
rfMidgetChrom = chromaticCenterWeight .* rfMidgetBipolar - ...
    midgetSurroundWeight .* chromaticSurroundWeight .* rfMidgetBipolar .* stageDelay(f,midgetSurroundDelay/1000);
rfMidgetLum = g.*rfMidgetLum;
rfMidgetChrom = g.*rfMidgetChrom;

end