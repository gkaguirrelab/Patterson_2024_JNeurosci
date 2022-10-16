%% fitRGCFResponseData
%
% Loads RGC temporal sensitivity data from Solomon et al (2002, 2005) and
% then fits these data in the complex fourier domain using a cascading
% low-pass filter model.
%
% The fitting is conducted simultaneously for parasol and midget responses
% at three eccentricities, and for the gain and phase components of the
% filter response.
%

%% Housekeeping
rng;
searchFlag = false;


%% Load the flicker response data
[midgetData, parasolData] = loadRGCResponseData();


%% Set model constants
eccFields = {'e0','e20','e30'};
% We assume that the RGC data were sampled from somewhere within these
% eccentricity ranges
eccBins = {[0 10],[20 30],[30 47]};
% The objective function attempts to fit both the gain and phase of the
% filter. The unwrapped phase units are large, so this is how they are
% scaled to be placed on a more equal footing with the gain.
phaseErrorScale = 1/400;
% The search includes a regularization that attempts to shrink the
% difference in some parameter values across eccentricity. This is the
% factor by which these differences are scaled in the objective function.
shrinkErrorScale = 20;


%% Set the p0 values
% There are 8 "block" parameters, and these 8 parameters are allowed to
% vary for each RGC class (midget and parasol) and each of the three
% eccentricities, leading to 6 blocks in total. The same p0 values are used
% to initialize all blocks
g = 4; % Overall gain
k = 0.67; % relative strength of the "lead compensators" (feedback stages)
cfInhibit = 20; % Corner frequency of the inhibitory stage
cf2ndStage = 40; % Corner frequency of 2nd order filter
Q = 1.0; % The "quality" parameter of the 2nd order filter
surroundWeight = 0.8667; % Weight of the surround relative to center
surroundDelay = 3; % Delay (in msecs) of the surround relative to center
eccProportion = 0.25; % Position within the eccentricity bin to calculate LM ratios

% There are 3 "common" parameters that define the temporal properties of
% the L and M cones, and the LM cone ratio.
cfCone = 15; % Corner frequency of the "cone" low-pass stage
coneDelay = 14; % Delay (in msecs) imposed at the "cone" stage
LMRatio = 1.0; % Ratio of L to M cones


%% Assemble the p0 and bounds
% We define here bounds on the block parameters
blockParamNames = {'g','k','cfInhibit', 'cf2ndStage', 'Q', 'surroundWeight', 'surroundDelay', 'eccProportion'};
p0Block = [g, k, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay, eccProportion];
lbBlock =  [2, 0.01, 01, 20, 0.25, 0.7, 01, 0.05];
plbBlock = [3, 0.40, 02, 30, 0.50, 0.8, 02, 0.10];
pubBlock = [5, 0.80, 40, 60, 2.50, 1.0, 05, 0.90];
ubBlock =  [6, 0.90, 60, 90, 3.00, 1.1, 10, 0.95];

% This vector controls for each of the block parameters if the shrink
% regularization is applied to match values across ecccentricity. A
% different set of choices are used for the midget and parasol models, as I
% find that the midget data can only be fit properly if greater flexibility
% across eccentricity is allowed.
shrinkParams = {...
    [false false false false false true false false], ... % midget shrink
    [true true true false false true true true] ...       % parasol shrink
    };

% Assemble the p0 and bounds out of the block and common parameters
p0 = [p0Block p0Block p0Block p0Block p0Block p0Block cfCone coneDelay LMRatio];
lb = [lbBlock lbBlock lbBlock lbBlock lbBlock lbBlock 10 10 0.50];
plb = [plbBlock plbBlock plbBlock plbBlock plbBlock plbBlock 12 12 0.90];
pub = [pubBlock pubBlock pubBlock pubBlock pubBlock pubBlock 18 15 1.10];
ub = [ubBlock ubBlock ubBlock ubBlock ubBlock ubBlock 20 18 2.00];

% Derive and set some values we will need later
nBlockParams = length(p0Block);
nEccBands = length(eccFields);
nCellClasses = 2;

% Here is a seed from a prior search with good performance
% fValGain: 4.53, fValPhase: 0.63, shrinkMidget: 0.00, shrinkParasol: 0.10 
p0 = [ 3.2832766771, 0.5869017720, 5.6224342804, 38.6179661751, 1.2765139937, 0.8945801139, 2.2711495161, 0.7119268656, 3.7975504398, 0.6717520714, 15.9822773454, 41.3446515799, 1.3080786467, 0.8946221888, 2.2795695066, 0.1550450563, 4.0466696024, 0.7272899985, 32.5966460044, 51.7450073361, 2.5190134645, 0.8946240008, 3.9494566917, 0.4883407354, 4.0075089931, 0.4785271406, 4.5361354048, 33.8026463985, 0.6534774899, 0.9341797292, 3.7295443416, 0.4988916636, 4.0461232066, 0.4786982656, 4.5361491743, 55.3475332260, 0.6540495753, 0.9343281448, 3.7301388979, 0.4995038033, 4.0255140662, 0.4789317966, 4.5378192406, 58.7253302336, 2.4650032520, 0.9344518602, 3.7305910289, 0.4999736309, 14.1357336044, 13.6896866262, 0.9967510045 ];


%% Define the objective and non-linear bound
myFit = @(p,verbose) rgcFitObjective(p,midgetData,parasolData,...
    shrinkParams,nBlockParams,eccFields,eccBins,...
    phaseErrorScale,shrinkErrorScale,verbose);
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


%% Report results
% Call the objective at the solution to report the fVals
myFit(p,true);

% Print the parameters in a format to be used as a seed in future searches
str = 'p0 = [ ';
for ss=1:length(p); str = [str sprintf('%2.10f, ',p(ss))]; end
str = [str(1:end-2) ' ];\n'];
fprintf(str);

% Reshape and store the parameters
LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:nBlockParams*nEccBands*nCellClasses),[nBlockParams,nEccBands,nCellClasses]);

% Report the common params
fprintf('cfCone: %2.2f, coneDelay: %2.2f, LMRatio: %2.2f \n',cfCone,coneDelay,LMRatio)


%% Interpolate params across eccentricity

% This is the position within each RGC recording bin that the model has
% identified as the correct location
eccDegs = zeros(1,nEccBands);
for ee = 1:nEccBands
    eccDegs(ee) = eccBins{ee}(1)+squeeze(p(8,ee,1)*range(eccBins{ee}));
end

% A simple linear interpolation and extrapolation, bounded by maximum
% returned values from the search
myInterpObj = @(v,xq,ii) max([repmat(plbBlock(ii),1,length(xq)); min([repmat(max(v),1,length(xq)); interp1(eccDegs,v,xq,'linear','extrap')])]);

% Loop across cells
for cc=1:2
    % Loop across the 7 params that vary with eccentricty
    for ii=1:nBlockParams-1
        % The values for this param across eccentricity, and the mean
        y = squeeze(p(ii,:,cc));
        pMean(ii,cc) = mean(y);
        % Fit the values and store the fit
        pFitByEccen{ii,cc} = @(xq) myInterpObj(y,xq,ii);
    end
end


%% Save the temporalModel
rgcTemporalModel.p = p;
rgcTemporalModel.LMRatio = LMRatio;
rgcTemporalModel.coneDelay = coneDelay;
rgcTemporalModel.cfCone = cfCone;
rgcTemporalModel.pMean = pMean;
rgcTemporalModel.pFitByEccen = pFitByEccen;
rgcTemporalModel.meta.blockParamNames = blockParamNames;
rgcTemporalModel.meta.eccFields = eccFields;
rgcTemporalModel.meta.eccBins = eccBins;
rgcTemporalModel.meta.lbBlock = lbBlock;
rgcTemporalModel.meta.ubBlock = ubBlock;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
save(savePath,'rgcTemporalModel');




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions


%% nonbcon
function c = nonbcon(p,nBlockParams,nEccBands)

% Sometimes BADS sends an empty set of parameters
if isempty(p)
    c=1;
end

% Enforce that some parameters, such as delay and filter frequency,
% increase in value across eccentricity
for ii=1:size(p,1)
    subP = p(ii,:);
    subP = reshape(subP(1:nBlockParams*nEccBands*2),[nBlockParams,nEccBands,2]);
    for cc=1:2
        tempP=squeeze(subP(:,:,cc));
        if ...
                any(diff(tempP(3,:))<0) || ... % force cfInhibit to increase with eccentricity
                any(diff(tempP(6,:))<0) || ... % force surroundWeight to increase with eccentricity
                any(diff(tempP(7,:))<0)        % force surroundDelay to increase with eccentricity
            c(ii,cc)=1;
        else
            c(ii,cc)=0;
        end
    end
end
c=any(c')';
end

%                 any(diff(tempP(4,:))<0) || ... % force cf2ndStage to increase with eccentricity
%                 any(diff(tempP(5,:))<0) || ... % force 2nd stage Q to increase with eccentricity


%% rgcFitObjective
function fVal = rgcFitObjective(p,midgetData,parasolData,shrinkParams,nBlockParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose)
% Given the midget and parasol data, as well as the full set of parameters
% across cell classes and eccentricity, derive the model fit error, with
% separate initial terms for fitting the gain and phase components of the
% transfer function. Additionally, we implement a "shrink" error that is
% used to drive certain sets of parameters to be the same across
% eccentricity locations.

% Extract and reshape the parameters
LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:nBlockParams*3*2),[nBlockParams,3,2]);

% Loop across eccentricity bands. Can make this parfor for speed
for ee = 1:length(eccFields)

    %% Midgets
    % Extract the parameter values for this eccentricity band
    pBlock = squeeze(p(:,ee,1));
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % The eccProportion parameter is used to determine the precise eccentricity
    % location within the range provided by the eccBin
    eccProportion = pBlock(8);
    eccDeg = eccBin(1)+eccProportion*(range(eccBin));

    % Get the temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum] = parseParamsMidget(pBlock, cfCone, coneDelay, LMRatio, eccDeg);

    % Derive the complex fourier domain TTF from the symbolic equations
    chromTTF = double(subs(rfMidgetChrom,midgetData.(eccField).chromatic.f));
    lumTTF = double(subs(rfMidgetLum,midgetData.(eccField).luminance.f));

    % Error in fitting the gain values
    chromGainErrorMidget(ee) = norm(midgetData.(eccField).chromatic.g - abs(chromTTF));
    lumGainErrorMidget(ee) = norm(midgetData.(eccField).luminance.g - abs(lumTTF));

    % Error in fitting the phase values
    chromPhaseErrorMidget(ee) = norm(midgetData.(eccField).chromatic.p - unwrap(angle(chromTTF))*(180/pi));
    lumPhaseErrorMidget(ee) = norm(midgetData.(eccField).luminance.p - unwrap(angle(lumTTF))*(180/pi));

    %% Parasols
    pBlock = squeeze(p(:,ee,2));

    % Get the temporal RF defined by these parameters
    rfParasaolLum = parseParamsParasol(pBlock, cfCone, coneDelay);

    % Derive the complex fourier domain TTF from the symbolic equations
    lumTTF = double(subs(rfParasaolLum,parasolData.(eccField).luminance.f));

    % Error in fitting the gain values    chromGainErrorMidget(ee) = norm(midgetData.(eccField).chromatic.g - abs(chromTTF));
    lumGainErrorParasol(ee) = norm(parasolData.(eccField).luminance.g - abs(lumTTF));

end

% Take the L2 norm of all of the model fits, scaling the phase values
fValGain = norm([lumGainErrorMidget chromGainErrorMidget lumGainErrorParasol]);
fValPhase = norm([lumPhaseErrorMidget chromPhaseErrorMidget])*phaseErrorScale;

% Calculate a regularization that attempts to match parameter values
% across eccentricity. Do this separately for the midget and parasol model
pSub = squeeze(p(shrinkParams{1},:,1));
fValShrinkMidget = shrinkErrorScale * norm(std(pSub,[],2)./mean(pSub,2));
pSub = squeeze(p(shrinkParams{2},:,2));
fValShrinkParasol = shrinkErrorScale * norm(std(pSub,[],2)./mean(pSub,2));

% Report the values
if verbose
    fprintf('fValGain: %2.2f, fValPhase: %2.2f, shrinkMidget: %2.2f, shrinkParasol: %2.2f \n',fValGain,fValPhase,fValShrinkMidget,fValShrinkParasol);
end

% Combine errors
fVal = fValGain + fValPhase + fValShrinkMidget + fValShrinkParasol;

% Catch and error out if there is a nan in here
if isnan(fVal)
    error('Encountered a nan value in the objective function');
end

end
