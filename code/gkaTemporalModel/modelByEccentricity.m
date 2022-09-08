clear

loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','temporalModel.mat');
load(loadPath,'temporalModel');

% Define a set of eccentricities in degrees, log spaced and centered within
% each of the bins of the Mt Sinai data
eccDegs = logspace(log10(1),log10(64),13); % eccentricities in degrees,

% A high-rez version of eD for plotting interpolated fits
eDFit = logspace(log10(1),log10(64),130);

% Drasdo 2007 equation for the midget fraction as a function of
% eccentricity
midgetFraction = @(r) 0.8928*(1+r/41.03).^(-1);


eccDegs = [20 25 30 35 40];

figHandleMidgetChrom = figure;
figHandleMidgetLum = figure;
figHandleParasolLum = figure;
figHandleMixLum = figure;

lineStyles = {'-','--','-.',':','-'};

for ee = 1:5

    pBlock = [];
    for ii = 1:7
        pBlock(ii) = temporalModel.pFitByEccen{ii,1}(eccDegs(ee));
    end

    [rfMidgetChrom, rfMidgetLum] = parseParamsMidget(pBlock, ...
        temporalModel.cfCone, temporalModel.coneDelay, temporalModel.LMRatio, eccDegs(ee));

    plotRF(rfMidgetChrom,figHandleMidgetChrom,[lineStyles{ee} 'r']);
    plotRF(rfMidgetLum,figHandleMidgetLum,[lineStyles{ee} 'k']);

    pBlock = [];
    for ii = 1:7
        pBlock(ii) = temporalModel.pFitByEccen{ii,2}(eccDegs(ee));
    end

    [rfParasolLum] = parseParamsParasol(pBlock, ...
        temporalModel.cfCone, temporalModel.coneDelay);

    plotRF(rfParasolLum,figHandleParasolLum,[lineStyles{ee} 'b']);

    rfMixLum{1} = midgetFraction(eccDegs(ee)).*rfMidgetLum;
    rfMixLum{2} = (1-midgetFraction(eccDegs(ee))).*rfParasolLum;
    plotRF(rfMixLum,figHandleMixLum,[lineStyles{ee} 'm']);


end

