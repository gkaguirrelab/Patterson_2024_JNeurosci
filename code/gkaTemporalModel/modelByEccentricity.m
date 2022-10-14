clear

loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','temporalModel.mat');
load(loadPath,'temporalModel');

% Define a set of eccentricities in degrees, log spaced and centered within
% each of the bins of the Mt Sinai data
eccDegBinEdges = logspace(log10(1),log10(90),20); % eccentricities in degrees,

% A high-rez version of eD for plotting interpolated fits
eDFit = logspace(log10(1),log10(64),130);

% Drasdo 2007 equation for the midget fraction as a function of
% eccentricity
midgetFraction = @(r) 0.8928*(1+r/41.03).^(-1);

figHandleMixLum = figure;
figHandleMixChrom = figure;

lineStyles = {'-','--','-.',':','.','o'};

rfMixChrom = {};
rfMixLum = {};

for ee = 1:length(eccDegBinEdges)-1

    eccDeg = mean([eccDegBinEdges(ee),eccDegBinEdges(ee+1)]);
    ringArea = pi*eccDegBinEdges(ee+1)^2 - pi*eccDegBinEdges(ee)^2;

    pBlock = [];
    for ii = 1:7
        pBlock(ii) = temporalModel.pFitByEccen{ii,1}(eccDeg);
    end

    [rfMidgetChrom, rfMidgetLum, f] = parseParamsMidget(pBlock, ...
        temporalModel.cfCone, temporalModel.coneDelay, temporalModel.LMRatio, eccDeg);

    rfMixChrom{end+1} = ringArea*midgetFraction(eccDeg).*rfMidgetChrom;

    pBlock = [];
    for ii = 1:7
        pBlock(ii) = temporalModel.pFitByEccen{ii,2}(eccDegBinEdges(ee));
    end

    [rfParasolLum] = parseParamsParasol(pBlock, ...
        temporalModel.cfCone, temporalModel.coneDelay);

    rfMixLum{end+1} = ringArea*midgetFraction(eccDeg).*rfMidgetLum;
    rfMixLum{end+1} = ringArea*(1-midgetFraction(eccDeg)).*rfParasolLum;


end


plotRF(rfMixChrom,figHandleMixChrom,'-r');
plotRF(rfMixLum,figHandleMixLum,'-k');
