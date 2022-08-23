% fitRGCFResponseData
%

rng;

%% Load the flicker response data
[midgetData, parasolData] = loadRGCResponseData();


%% Set model constants
eccFields = {'e0','e20','e30'};
eccBins = {[0 10],[20 30],[30 47]};
phaseErrorScale = 1/400;
shrinkErrorScale = 20;


%% Set the p0 values
g = 4; % Overall gain
k = 0.67; % relative strength of the "lead compensators" (feedback stages)
cfLowPass = 20; % Corner frequency of the "bipolar" low-pass stage
cfInhibit = 20; % Corner frequency of the inhibitory stage
cf2ndStage = 40; % Corner frequency of 2nd order filter
Q = 1.0; % The "quality" parameter of the 2nd order filter
surroundWeight = 0.8667; % Weight of the surround relative to center
surroundDelay = 3; % Delay (in msecs) of the surround relative to center
eccProportion = 0.25; % Position within the eccentricity bin to calculate LM ratios

cfCone = 12; % Corner frequency of the "cone" low-pass stage
coneDelay = 14; % Delay (in msecs) impposed by the "cone" stage
LMRatio = 1.0; % Ratio of L to M cones


%% Define p0 and bounds
p0Block = [g, k, cfLowPass, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay, eccProportion];
lbBlock =  [0, 0.40, 05, 05, 020, 0.50, 0.0, 01, 0.05];
plbBlock = [3, 0.50, 10, 10, 050, 0.75, 0.5, 02, 0.10];
pubBlock = [6, 0.85, 40, 70, 090, 2.50, 1.0, 08, 0.90];
ubBlock = [10, 0.90, 50, 80, 100, 3.00, 1.5, 10, 0.95];
shrinkParams = [false false true false false false true false false];

p0 = [p0Block p0Block p0Block cfCone coneDelay LMRatio];
lb = [lbBlock lbBlock lbBlock 05 5 0.1];
plb = [plbBlock plbBlock plbBlock 10 10 0.33];
pub = [pubBlock pubBlock pubBlock 20 20 3];
ub = [ubBlock ubBlock ubBlock 25 25 10];

% Could replace the default p0 here with a seed from a prior search
p0 = [ 3.2428, 0.5920, 20.2840, 6.5189, 41.4429, 1.1746, 0.8646, 2.1941, 0.6054, 3.8758, 0.6809, 20.2869, 18.3447, 41.4600, 1.1795, 0.8647, 2.1946, 0.0501, 4.0142, 0.7000, 20.2882, 31.0958, 51.7202, 2.4805, 0.8647, 4.3085, 0.2639, 11.3229, 13.7054, 1.0001 ];


%% Define the objective
myFit = @(p,verbose) rgcFitObjective(p,midgetData,shrinkParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose);
myObj = @(p) myFit(p,false);


%% Options
% Our objective function is deterministic
options.UncertaintyHandling = 0;


%% Search
p = bads(myObj,p0,lb,ub,plb,pub,@nonbcon,options); 

% Call the objective at the solution to report the fVals
myFit(p,true);

% Print the parameters in a format to be used as a seed in future searches 
str = 'p0 = [ ';
for ss=1:length(p); str = [str sprintf('%2.4f, ',p(ss))]; end
str = [str(1:end-2) ' ];\n'];
fprintf(str);


%% Report the results
LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:27),[9,3]);

% Report the common params
fprintf('cfCone: %2.2f, coneDelay: %2.2f, LMRatio: %2.2f \n',cfCone,coneDelay,LMRatio)

% Dump out the reshaped p values
p

% Plot the temporal receptive fields at the "cone" and "bipolar" stages
figHandle = figure();


% Plot each eccentricity band
for ee = 1:3

    % Extract the parameter values for this eccentricity band
    pBlock = p(:,ee);
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % Get the temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum, eccDegs, rfLMCone, rfBipolar] = ...
        parseParams(pBlock, LMRatio, cfCone, coneDelay, eccBin);

    % Plot the temporal RFs
    figHandle = figure();
    plotRF(rfMidgetLum,figHandle,'-k');
    plotRF(rfMidgetChrom,figHandle,'-r');
    plotRF(rfLMCone,figHandle,'-g',3);
    subplot(3,1,1);
    loglog(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.g,'*r');
    loglog(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.g,'*k');
    title(sprintf('Eccentricity = %2.1f',eccDegs));
    subplot(3,1,2);
    semilogx(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.p,'*r');
    semilogx(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.p,'*k');

end


% Local functions

function c = nonbcon(p)

if isempty(p)
    c=1;
end

% Enforce that some parameters, such as delay and filter frequency,
% increase or decrease in value across eccentricity
for ii=1:size(p,1)
    tempP = reshape( squeeze(p(ii,1:27)),[9,3]);
    if ...
            any(diff(tempP(3,:))<0) || ... % force cfLowPass to increase with eccentricity
            any(diff(tempP(4,:))<0) || ... % force cfInhibit to increase with eccentricity
            any(diff(tempP(5,:))<0) || ... % force cf2ndStage to increase with eccentricity
            any(diff(tempP(6,:))<0) || ... % force 2nd stage Q to increase with eccentricity
            any(diff(tempP(7,:))<0) || ... % force surroundWeight to increase with eccentricity
            any(diff(tempP(8,:))<0)        % force surroundDelay to increase with eccentricity
        c(ii)=1;
    else
        c(ii)=0;
    end
end
c=c';
end



function [rfMidgetChrom, rfMidgetLum, eccDegs, rfLMCone, rfBipolar] = parseParams(pBlock, LMRatio, cfCone, coneDelay, eccBin)

g = pBlock(1); k = pBlock(2);
cfLowPass = pBlock(3); cfInhibit = pBlock(4);
cf2ndStage = pBlock(5); Q = pBlock(6);
surroundWeight = pBlock(7); surroundDelay = pBlock(8);
eccProportion = pBlock(9);

eccDegs = eccBin(1)+eccProportion*(range(eccBin));

tmpCenterWeight = [];
tmpSurroundWeight = [];

for cc = 1:1000
    [tmpC,tmpS] = ...
        returnChromaticWeights(eccDegs,LMRatio);
    tmpCenterWeight(cc) = tmpC;
    tmpSurroundWeight(cc) = tmpS;
end
chromaticCenterWeight = mean(abs(tmpCenterWeight));
chromaticSurroundWeight = mean(abs(tmpSurroundWeight));

[rfLMCone, rfBipolar, rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
    cfCone, coneDelay, ...
    g, k, cfLowPass, cfInhibit, cf2ndStage, Q, ...
    surroundWeight, surroundDelay, ...
    chromaticCenterWeight, chromaticSurroundWeight);

end


function fVal = rgcFitObjective(p,midgetData,shrinkParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose)

LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:27),[9,3]);

chromGainError = [];
lumGainError = [];
chromPhaseError = [];
lumPhaseError = [];

% Loop across eccentricity bands
parfor ee = 1:3

    % Extract the parameter values for this eccentricity band
    pBlock = p(:,ee);
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % Get the temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum] = parseParams(pBlock, LMRatio, cfCone, coneDelay, eccBin);

    % Derive the complex fourier domain TTF from the symbolic equations
    chromTTF = double(subs(rfMidgetChrom,midgetData.(eccField).chromatic.f));
    lumTTF = double(subs(rfMidgetLum,midgetData.(eccField).luminance.f));

    % Error in fitting the gain values
    chromGainError(ee) = norm(midgetData.(eccField).chromatic.g - abs(chromTTF));
    lumGainError(ee) = norm(midgetData.(eccField).luminance.g - abs(lumTTF));

    % Error in fitting the phase values
    chromPhaseError(ee) = norm(midgetData.(eccField).chromatic.p - unwrap(angle(chromTTF))*(180/pi));
    lumPhaseError(ee) = norm(midgetData.(eccField).luminance.p - unwrap(angle(lumTTF))*(180/pi));

end

% Take the L2 norm of all of the model fits, scaling the phase values
fValGain = norm([lumGainError chromGainError]);
fValPhase = norm([lumPhaseError chromPhaseError])*phaseErrorScale;

% Calculate a regularization that attempts to match parameter values
% across eccentricity
fValShrink = shrinkErrorScale * norm(std(p(shrinkParams,:),[],2)./mean(p(shrinkParams,:),2));

% Report the values
if verbose
    fprintf('fValGain: %2.2f, fValPhase: %2.2f, shrink: %2.2f \n',fValGain,fValPhase,fValShrink);
end

% Combine errors
fVal = fValGain + fValPhase + fValShrink;

if isnan(fVal)
    error('Encountered a nan value in the objective function');
end

end
