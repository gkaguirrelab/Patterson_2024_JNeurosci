loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','midgetModel.mat');
load(loadPath,'midgetModel');

% Define a set of eccentricities in degrees, log spaced and centered within
% each of the bins of the Mt Sinai data
eD = logspace(log10(1),log10(64),13); % eccentricities in degrees,

% A high-rez version of eD for plotting interpolated fits
eDFit = logspace(log10(1),log10(64),130);


figHandleChrom = figure;
figHandleLum = figure;

lineStyles = {'-','--','-.',':'};
eccVals = [6.25 12.5 25 50];
for ee = 1:4
[rfMidgetChrom, rfMidgetLum] = parseParams(midgetModel,eccVals(ee));
plotRF(rfMidgetChrom,figHandleChrom,[lineStyles{ee} 'r']);
plotRF(rfMidgetLum,figHandleLum,[lineStyles{ee} 'k']);
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
function [rfMidgetChrom, rfMidgetLum] = parseParams(midgetModel,eccDegs)

% Derive the center and surround chromatic weights
tmpCenterWeight = [];
tmpSurroundWeight = [];

LMRatio = midgetModel.LMRatio;

for cc = 1:1000
    [tmpC,tmpS] = ...
        returnChromaticWeights(eccDegs,LMRatio);
    tmpCenterWeight(cc) = tmpC;
    tmpSurroundWeight(cc) = tmpS;
end
chromaticCenterWeight = mean(abs(tmpCenterWeight));
chromaticSurroundWeight = mean(abs(tmpSurroundWeight));

% Obtain the model parameters for this eccentricity and assign to a
% variable in the current function
sigmoidFit = @(x) x(1) + x(2) - x(2)*exp( - (eccDegs./x(3)).^x(4) );
for ii=1:7
    p(ii) = sigmoidFit(midgetModel.pByEccExpParams(ii,:));
    assignParam(midgetModel.blockParamNames{ii},p(ii));
end

Q = min([Q 2.5076]);
cf2ndStage = min([cf2ndStage 52.083]);
surroundDelay = min([surroundDelay 4.36621]);

% Create the RFs
[~, ~, rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
    midgetModel.cfCone, midgetModel.coneDelay, ...
    g, k, cfInhibit, cf2ndStage, Q, ...
    surroundWeight, surroundDelay, ...
    chromaticCenterWeight, chromaticSurroundWeight);

end


function assignParam(name,value)
    assignin('caller',name,value);
end