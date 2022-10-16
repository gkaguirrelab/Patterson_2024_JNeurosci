function [rfPostRetinalMidgetChrom,rfPostRetinalLum,...
    rfPostRetinalMidgetLum,rfPostRetinalParasolLum,...
    rfRGCMidgetChrom,rfRGCMidgetLum,rfRGCParasolLum,...
    rfMidgetBipolar,rfParasolBipolar,rfLMCone] = ...
        returnPostRetinalResponses(eccDeg,secondOrderFc,secondOrderQ,surroundIndex,surroundDelay,nSubtractions)

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Drasdo 2007 equation for the midget fraction as a function of
% eccentricity
midgetFraction = @(eccDeg) 0.8928*(1+eccDeg/41.03).^(-1);

% Loop over the meridians and obtain the RGCf density functions
meridianAngles = [0 90 180 270];
for mm = 1:4
    totalRGCfDensity{mm} = watsonTotalRFDensityByEccenDegVisual(meridianAngles(mm));
end

% Obtain the mean, total ganglion receptive field density as a function of
% eccentricity
totalRGCfDensityAtEcc = @(eccDeg) mean(cellfun(@(x) x(eccDeg),totalRGCfDensity));

% Get the rfs for midget and parasol luminance at this eccentricity
[rfPostRetinalMidgetChrom,rfPostRetinalMidgetLum,rfPostRetinalParasolLum,...
    rfRGCMidgetChrom,rfRGCMidgetLum,rfRGCParasolLum,...
    rfMidgetBipolar,rfParasolBipolar,rfLMCone] = ...
    returnCellModel(rgcTemporalModel,...
    eccDeg,secondOrderFc,secondOrderQ,...
    surroundIndex,surroundDelay,nSubtractions);

% Scale by the density of RGCf at this eccentricity, accounting for the
% change in surface area as a function of eccentricity (eccDeg/90), and
% then scale by the midget fraction.
rfPostRetinalMidgetChrom = rfPostRetinalMidgetChrom ...
    .* totalRGCfDensityAtEcc(eccDeg) .* (eccDeg/90) ...
    .* midgetFraction(eccDeg);

rfPostRetinalMidgetLum = rfPostRetinalMidgetLum ...
    .* totalRGCfDensityAtEcc(eccDeg) .* (eccDeg/90) ...
    .* midgetFraction(eccDeg);

rfPostRetinalParasolLum = rfPostRetinalParasolLum ...
    .* totalRGCfDensityAtEcc(eccDeg) .* (eccDeg/90) ...
    .* (1-midgetFraction(eccDeg));

% Combine the two luminance components in a single model
rfPostRetinalLum = {rfPostRetinalMidgetLum,rfPostRetinalParasolLum};

end



%% LOCAL FUNCTIONS
function [rfPostRetinalMidgetChrom,rfPostRetinalMidgetLum,rfPostRetinalParasolLum,...
    rfRGCMidgetChrom,rfRGCMidgetLum,rfRGCParasolLum,...
    rfMidgetBipolar,rfParasolBipolar,rfLMCone] = ...
        returnCellModel(temporalModel,eccDeg,secondOrderFc,secondOrderQ,surroundIndex,surroundDelay,nSubtractions)

syms f

% Obtain the RGC model parameters for this eccentricity
pBlockMidget = [];
pBlockParasol = [];
for ii = 1:7
    pBlockMidget(ii) = temporalModel.pFitByEccen{ii,1}(eccDeg);
    pBlockParasol(ii) = temporalModel.pFitByEccen{ii,2}(eccDeg);
end

% Derive the RGC models
[rfRGCMidgetChrom, rfRGCMidgetLum, ~, rfMidgetBipolar] = parseParamsMidget(pBlockMidget, ...
    temporalModel.cfCone, temporalModel.coneDelay, temporalModel.LMRatio, eccDeg);

[rfRGCParasolLum, rfLMCone, rfParasolBipolar] = parseParamsParasol(pBlockParasol, ...
    temporalModel.cfCone, temporalModel.coneDelay);

% Copy the RGC model into the post-retinal variables
rfPostRetinalMidgetChrom = rfRGCMidgetChrom; rfPostRetinalMidgetLum = rfRGCMidgetLum;
rfPostRetinalParasolLum = rfRGCParasolLum;

% Apply iterative, delayed surround subtraction. The expectation is one
% iteration for LGN, two iterations for V1
for ii = 1:nSubtractions
    rfPostRetinalMidgetChrom = rfPostRetinalMidgetChrom - surroundIndex*rfPostRetinalMidgetChrom.*stageDelay(f,surroundDelay/1000);
    rfPostRetinalMidgetLum = rfPostRetinalMidgetLum - surroundIndex*rfPostRetinalMidgetLum.*stageDelay(f,surroundDelay/1000);
    rfPostRetinalParasolLum = rfPostRetinalParasolLum - surroundIndex*rfPostRetinalParasolLum.*stageDelay(f,surroundDelay/1000);
end

% Apply a post-lgn, second order low pass filter
if nSubtractions == 2
    rfPostRetinalMidgetChrom = rfPostRetinalMidgetChrom.*stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
    rfPostRetinalMidgetLum = rfPostRetinalMidgetLum.*stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
    rfPostRetinalParasolLum = rfPostRetinalParasolLum.*stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
end

end