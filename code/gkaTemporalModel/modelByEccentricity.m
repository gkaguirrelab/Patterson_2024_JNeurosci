loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','temporalModel.mat');
load(loadPath,'temporalModel');

% Define a set of eccentricities in degrees, log spaced and centered within
% each of the bins of the Mt Sinai data
eD = logspace(log10(1),log10(64),13); % eccentricities in degrees,

% A high-rez version of eD for plotting interpolated fits
eDFit = logspace(log10(1),log10(64),130);

eccDegs = [5 10 20 40];
LMRatios = 1./[1 2 4 8];
figHandleChrom = figure;
%figHandleLum = figure;

lineStyles = {'-','--','-.',':'};
eccVals = [6.25 12.5 25 50];
for ee = 1:4

    pBlock = temporalModel.pMean(:,1);

    [rfMidgetChrom, rfMidgetLum] = parseParamsMidget(pBlock, ...
        temporalModel.cfCone, temporalModel.coneDelay, LMRatios(ee), eccDegs(1));
        
    plotRF(rfMidgetChrom,figHandleChrom,[lineStyles{ee} 'r']);
    %plotRF(rfMidgetLum,figHandleLum,[lineStyles{ee} 'k']);

end

