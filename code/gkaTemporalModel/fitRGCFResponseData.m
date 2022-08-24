% fitRGCFResponseData
%

rng;
searchFlag = false;


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
cfInhibit = 20; % Corner frequency of the inhibitory stage
cf2ndStage = 40; % Corner frequency of 2nd order filter
Q = 1.0; % The "quality" parameter of the 2nd order filter
surroundWeight = 0.8667; % Weight of the surround relative to center
surroundDelay = 3; % Delay (in msecs) of the surround relative to center
eccProportion = 0.25; % Position within the eccentricity bin to calculate LM ratios

cfCone = 15; % Corner frequency of the "cone" low-pass stage
coneDelay = 14; % Delay (in msecs) impposed by the "cone" stage
LMRatio = 1.0; % Ratio of L to M cones


%% Define p0 and bounds
blockParamNames = {'g','k','cfInhibit', 'cf2ndStage', 'Q', 'surroundWeight', 'surroundDelay', 'eccProportion'};
p0Block = [g, k, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay, eccProportion];
lbBlock =  [0, 0.40, 05, 020, 0.50, 0.0, 01, 0.05];
plbBlock = [3, 0.50, 10, 050, 0.75, 0.5, 02, 0.10];
pubBlock = [6, 0.85, 70, 090, 2.50, 1.0, 08, 0.90];
ubBlock = [10, 0.90, 80, 100, 3.00, 1.5, 10, 0.95];
shrinkParams = [false false false false false true false false];

p0 = [p0Block p0Block p0Block cfCone coneDelay LMRatio];
lb = [lbBlock lbBlock lbBlock 05 5 0.1];
plb = [plbBlock plbBlock plbBlock 10 10 0.33];
pub = [pubBlock pubBlock pubBlock 20 20 3];
ub = [ubBlock ubBlock ubBlock 25 25 10];

nBlockParams = length(p0Block);
nEccBands = length(eccFields);

% Here is a seed from a prior search with good performance
% fValGain: 1.60, fValPhase: 0.67, shrink: 0.00
p0 = [ 3.4397165179, 0.5859245405, 6.7901217937, 40.0400388241, 1.1437500417, 0.8688447177, 2.0883977413, 0.7998218775, 3.8780929148, 0.6806570783, 19.3246120214, 40.4715168476, 1.2194156796, 0.8688532561, 2.1265215874, 0.0543420553, 4.0079992712, 0.7113094196, 31.2511581182, 51.9151723385, 2.5076012462, 0.8688690066, 4.3271616101, 0.3479259253, 15.1120990515, 13.7365531921, 1.0085806677 ];


%% Define the objective and non-linear bound
myFit = @(p,verbose) rgcFitObjective(p,midgetData,shrinkParams,nBlockParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose);
myObj = @(p) myFit(p,false);
myNonbcon = @(p) nonbcon(p,nBlockParams,nEccBands);


%% Options
% The objective function is deterministic
optionsBADS.UncertaintyHandling = 0;


%% Search
if searchFlag
    p = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);
else
    p = p0;
end

% Call the objective at the solution to report the fVals
myFit(p,true);

% Print the parameters in a format to be used as a seed in future searches
str = 'p0 = [ ';
for ss=1:length(p); str = [str sprintf('%2.10f, ',p(ss))]; end
str = [str(1:end-2) ' ];\n'];
fprintf(str);


%% Report the results
LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:nBlockParams*3),[nBlockParams,nEccBands]);

% Report the common params
fprintf('cfCone: %2.2f, coneDelay: %2.2f, LMRatio: %2.2f \n',cfCone,coneDelay,LMRatio)

% Dump out the reshaped p values
p

% Plot each eccentricity band
for ee = 1:nEccBands

    % Extract the parameter values for this eccentricity band
    pBlock = p(:,ee);
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % Get the temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum, eccDegs(ee), rfLMCone, rfBipolar] = ...
        parseParams(pBlock, LMRatio, cfCone, coneDelay, eccBin);

    % Plot the temporal RFs
    figHandle = figure();
    plotRF(rfMidgetLum,figHandle,'-k');
    plotRF(rfMidgetChrom,figHandle,'-r');
    plotRF(rfLMCone,figHandle,'-g',3);
    subplot(3,1,1);
    loglog(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.g,'*r');
    loglog(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.g,'*k');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    subplot(3,1,2);
    semilogx(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.p,'*r');
    semilogx(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.p,'*k');

end

% Plot the parameters vs. eccentricity and obtain params x eccentricity
figure
pByEccExpParams = nan(nBlockParams,4);
optionsFMINCON = optimoptions('fmincon','Display','off');

sigmoidFit = @(x,support) x(1) + x(2) - x(2)*exp( - (support./x(3)).^x(4) );

for ii=1:nBlockParams-1
    subplot(4,2,ii)
    y = p(ii,:);
    plot(eccDegs,y,'ok')
    hold on
    x0 = [y(1),range(y),mean(eccDegs),10];
    myPlotFitObj = @(x) norm(sigmoidFit(x,eccDegs)-y);
    pByEccExpParams(ii,:)=fmincon(myPlotFitObj,x0,[],[],[],[],[],[],[],optionsFMINCON);
    eccDegsFit = 1:max(eccDegs)*1.05;
    plot(eccDegsFit,sigmoidFit(pByEccExpParams(ii,:),eccDegsFit),'-r')
    title(blockParamNames{ii})
    xlabel('Eccentricity [deg]'); ylabel('param value')
    if ii==paramToHold
        ylim([0 1]);
    end
end


%% Save the midgetModel structure
midgetModel.pByEccExpParams = pByEccExpParams;
midgetModel.blockParamNames = blockParamNames;
midgetModel.LMRatio = LMRatio;
midgetModel.coneDelay = coneDelay;
midgetModel.cfCone = cfCone;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','midgetModel.mat');
save(savePath,'midgetModel');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

function c = nonbcon(p,nBlockParams,nEccBands)

if isempty(p)
    c=1;
end

% Enforce that some parameters, such as delay and filter frequency,
% increase or decrease in value across eccentricity
for ii=1:size(p,1)
    tempP = reshape( squeeze(p(ii,1:nBlockParams*nEccBands)),[nBlockParams,nEccBands]);
    if ...
            any(diff(tempP(3,:))<0) || ... % force cfInhibit to increase with eccentricity
            any(diff(tempP(4,:))<0) || ... % force cf2ndStage to increase with eccentricity
            any(diff(tempP(5,:))<0) || ... % force 2nd stage Q to increase with eccentricity
            any(diff(tempP(6,:))<0) || ... % force surroundWeight to increase with eccentricity
            any(diff(tempP(7,:))<0)        % force surroundDelay to increase with eccentricity
        c(ii)=1;
    else
        c(ii)=0;
    end
end
c=c';
end



function [rfMidgetChrom, rfMidgetLum, eccDegs, rfLMCone, rfBipolar] = parseParams(pBlock, LMRatio, cfCone, coneDelay, eccBin)

g = pBlock(1); k = pBlock(2);
cfInhibit = pBlock(3); cf2ndStage = pBlock(4); Q = pBlock(5);
surroundWeight = pBlock(6); surroundDelay = pBlock(7);
eccProportion = pBlock(8);

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
    g, k, cfInhibit, cf2ndStage, Q, ...
    surroundWeight, surroundDelay, ...
    chromaticCenterWeight, chromaticSurroundWeight);

end


function fVal = rgcFitObjective(p,midgetData,shrinkParams,nBlockParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose)


LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:nBlockParams*3),[nBlockParams,3]);

chromGainError = [];
lumGainError = [];
chromPhaseError = [];
lumPhaseError = [];

% Loop across eccentricity bands
parfor ee = 1:length(eccFields)

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
