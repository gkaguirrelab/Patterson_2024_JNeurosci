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
searchFlag = true;


%% Load the flicker response data
rgcData = loadRGCResponseData();


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
lbBlock =  [0.00, 0.01, 01, 20, 0.25, 0.7, 01, 0.05];
plbBlock = [0.10, 0.40, 02, 30, 0.50, 0.8, 02, 0.10];
pubBlock = [1.25, 0.80, 40, 60, 2.50, 1.0, 05, 0.90];
ubBlock =  [2.00, 0.90, 60, 90, 3.00, 1.1, 10, 0.95];

% This vector controls for each of the block parameters if the shrink
% regularization is applied to match values across ecccentricity. A
% different set of choices are used for the midget and parasol models, as I
% find that the midget data can only be fit properly if greater flexibility
% across eccentricity is allowed.
shrinkParams = {...
    [false false false false false true false false], ... % midget shrink
    [true true true false false true true true] ...       % parasol shrink
    [true true true true true true true true] ...         % bistratified shrink
    };

% Assemble the p0 and bounds out of the block and common parameters
nBlocks = 6;
p0 = [repmat(p0Block,1,nBlocks) cfCone coneDelay LMRatio];
lb = [repmat(lbBlock,1,nBlocks) 10 10 0.50];
plb = [repmat(plbBlock,1,nBlocks) 12 12 0.90];
pub = [repmat(pubBlock,1,nBlocks) 18 15 1.10];
ub = [repmat(ubBlock,1,nBlocks) 20 18 2.00];

% Derive and set some values we will need later
nBlockParams = length(p0Block);
nEccBands = length(eccFields);
nCellClasses = 2;

% Here is a seed from a prior search with good performance
% fValGain: 4.55, fValPhase: 0.66, shrinkMidget: 0.00, shrinkParasol: 0.03 
p0 = [ ...
    0.1972097355, 0.5846308186, 5.0149334709, 39.8756076121, 1.5605660335, 0.8973745150, 2.4044085092, 0.7469553411, ...
    0.5832596217, 0.6475761358, 16.3549848095, 42.5415901253, 1.3650369485, 0.8973755821, 2.4506626798, 0.0500003678, ...
    1.1184292023, 0.7163389747, 32.4052113911, 51.3445703584, 2.5330080011, 0.8974094916, 3.8386941140, 0.7034826188, ...
    1.0643874637, 0.4814463600, 4.4342387618, 34.3600976944, 0.6029321683, 0.9352258071, 3.7132009941, 0.5003314620, ...
    1.0654956466, 0.4808840309, 4.4345761549, 58.0414818227, 0.7000975664, 0.9353066738, 3.7183583268, 0.5005013968, ...
    1.0647037918, 0.4813159307, 4.4349439070, 59.1540781916, 2.4885065617, 0.9353385903, 3.7192254514, 0.5011164645, ...
    14.1423620647, 13.6251217899, 0.9929001757 ];

%% Define the objective and non-linear bound
myFit = @(p,verbose) rgcFitObjective(p,rgcData,...
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
for ss=1:length(p)
    str = [str sprintf('%2.10f, ',p(ss))];
    if mod(ss,nBlockParams)==0
        str = [str '...\n'];
    end
end
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
function fVal = rgcFitObjective(p,rgcData,shrinkParams,nBlockParams,eccFields,eccBins,phaseErrorScale,shrinkErrorScale,verbose)
% Given the midget and parasol data, as well as the full set of parameters
% across cell classes and eccentricity, derive the model fit error, with
% separate initial terms for fitting the gain and phase components of the
% transfer function. Additionally, we implement a "shrink" error that is
% used to drive certain sets of parameters to be the same across
% eccentricity locations.

nEccBins = length(eccFields);
nCellClasses = 2;

% Extract and reshape the parameters
LMRatio = p(end);
coneDelay = p(end-1);
cfCone = p(end-2);
p = reshape(p(1:nBlockParams*3*2),[nBlockParams,nEccBins,nCellClasses]);

% Loop across eccentricity bands. Can make this parfor for speed
for ee = 1:length(eccFields)

    % The eccentricity domain
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    %% Midget
    % Extract the parameter values for this eccentricity band
    pBlock = squeeze(p(:,ee,1));

    % The eccProportion parameter is used to determine the precise
    % eccentricity location within the range provided by the eccBin. This
    % is needed for the midget calculation (as it is sensitivity to the
    % random cone wiring at this position) but not the other RGC classe.
    eccProportion = pBlock(8);
    eccDeg = eccBin(1)+eccProportion*(range(eccBin));

    % Get the temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum] = parseParamsMidget(pBlock, cfCone, coneDelay, LMRatio, eccDeg);

    % Derive the complex fourier domain TTF from the symbolic equations
    chromTTF = double(subs(rfMidgetChrom,rgcData.midget.(eccField).chromatic.f));
    lumTTF = double(subs(rfMidgetLum,rgcData.midget.(eccField).luminance.f));

    % Error in fitting the gain values
    chromGainErrorMidget(ee) = norm(rgcData.midget.(eccField).chromatic.g - abs(chromTTF));
    lumGainErrorMidget(ee) = norm(rgcData.midget.(eccField).luminance.g - abs(lumTTF));

    % Error in fitting the phase values
    chromPhaseErrorMidget(ee) = norm(rgcData.midget.(eccField).chromatic.p - unwrap(angle(chromTTF))*(180/pi));
    lumPhaseErrorMidget(ee) = norm(rgcData.midget.(eccField).luminance.p - unwrap(angle(lumTTF))*(180/pi));

    %% Parasol
    % Extract the parameter values for this eccentricity band
    pBlock = squeeze(p(:,ee,2));

    % Get the temporal RF defined by these parameters
    rfParasaolLum = parseParamsParasol(pBlock, cfCone, coneDelay);

    % Derive the complex fourier domain TTF from the symbolic equations
    lumTTF = double(subs(rfParasaolLum,rgcData.parasol.(eccField).luminance.f));

    % Error in fitting the gain values    chromGainErrorMidget(ee) = norm(rgcData.midget.(eccField).chromatic.g - abs(chromTTF));
    lumGainErrorParasol(ee) = norm(rgcData.parasol.(eccField).luminance.g - abs(lumTTF));

    %% Bistratified
    % Extract the parameter values for this eccentricity band
    

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
