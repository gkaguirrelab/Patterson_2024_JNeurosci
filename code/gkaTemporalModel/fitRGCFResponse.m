function rgcTemporalModel = fitRGCFResponse(rgcSearchFlag,verboseFlag,paramSearch)

% Loads RGC temporal sensitivity data from Solomon et al (2002, 2005) and
% then fits these data in the complex fourier domain using a cascading
% low-pass filter model.
%
% The fitting is conducted simultaneously for parasol, midget, and
% bistratified responses at three eccentricities, and for the gain and
% phase components of the filter response.
%
% Examples:
%{
    % Create and save the temporal model results
    fitRGCFResponse(false,true);
%}

%% Housekeeping
rng;
if nargin==0
    rgcSearchFlag = false;
    verboseFlag = false;
    paramSearch = 'full';
end
if nargin==1
    verboseFlag = false;
    paramSearch = 'full';
end
if nargin==2
    paramSearch = 'full';
end


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

% This is the order in which the block parameters are stored across cell
% classes
cellClassIndices = {'midget','parasol','bistratified'};

% This vector controls for each of the block parameters if the shrink
% regularization is applied to match values across ecccentricity. A
% different set of choices are used for the midget and parasol models, as I
% find that the midget data can only be fit properly if greater flexibility
% across eccentricity is allowed.
shrinkParams = {...
    [false false false false false true false false], ... % midget shrink
    [true true true false false true true true] ...       % parasol shrink
    [true true true false false true true true] ...       % bistratified shrink
    };

% Derive and set some values we will need later
nBlockParams = length(blockParamNames);
nEccBands = length(eccFields);
nCellClasses = length(fieldnames(rgcData));

% Assemble the p0 and bounds out of the block and common parameters
nBlocks = nEccBands*nCellClasses;
p0 = [repmat(p0Block,1,nBlocks) cfCone coneDelay LMRatio];
lb = [repmat(lbBlock,1,nBlocks) 10 10 0.50];
plb = [repmat(plbBlock,1,nBlocks) 12 12 0.90];
pub = [repmat(pubBlock,1,nBlocks) 18 15 1.10];
ub = [repmat(ubBlock,1,nBlocks) 20 18 2.00];

% Here is a seed from a prior search with good performance
% midget - fValGain: 1.56, fValPhase: 0.66, fValShrink: 0.00
% parasol - fValGain: 3.25, fValPhase: 0.00, fValShrink: 2.58
% bistratified - fValGain: 0.46, fValPhase: 0.62, fValShrink: 0.00
p0 = [ ...
    0.2469418824, 0.5922375411, 5.8056983257, 39.0848216414, 1.3149221987, 0.8983960792, 2.1111500710, 0.7764686763, ...
    0.6415928598, 0.6496570826, 17.3118880045, 42.4663314968, 1.3954235017, 0.8984520018, 2.1172960848, 0.0508179903, ...
    1.0906113461, 0.7230467349, 30.8115966851, 51.7642246932, 2.5439632535, 0.8984555170, 3.7704597861, 0.7612058520, ...
    0.8828694341, 0.8999820508, 2.1840414273, 27.4735044287, 1.1172341323, 0.9387978947, 3.2146946587, 0.8726685094, ...
    1.0597182026, 0.8947987072, 2.1840659747, 55.5393972857, 0.7654854182, 0.9387979018, 3.5752151019, 0.9474541004, ...
    1.0453085768, 0.8447868041, 2.1840736126, 58.0722041866, 2.6753874166, 0.9387979452, 3.5752221443, 0.9412647515, ...
    0.8968017565, 0.4126953134, 9.1694768937, 20.9763595134, 0.4257805828, 0.9360351530, 3.1059570000, 0.0625000003, ...
    0.8968017565, 0.4126953134, 9.1694768937, 20.9763595134, 0.4257805828, 0.9360351530, 3.1059570000, 0.0625000003, ...
    0.8968017565, 0.4126953134, 9.1694768937, 20.9763595134, 0.4257805828, 0.9360351530, 3.1059570000, 0.0625000003, ...
    14.2376403809, 13.7579040527, 0.9915176392 ];


% The search is unable to operate across all the parameters at once. This
% code allows us to lock parameters and search for a single RGC class at a
% time
totalBlockParams = nBlockParams*nEccBands*nCellClasses;
switch paramSearch
    case 'midget'
        lockIdx = [nBlockParams*nEccBands+1:totalBlockParams, totalBlockParams+1:totalBlockParams+3];
    otherwise
        lockIdx = [];
end
lb(lockIdx) = p0(lockIdx); plb(lockIdx) = p0(lockIdx);
ub(lockIdx) = p0(lockIdx); pub(lockIdx) = p0(lockIdx);


%% Define the objective and non-linear bound
myFit = @(pRGC,verbose) rgcFitOverallObjective(pRGC,rgcData,...
    eccFields,eccBins,nBlockParams,cellClassIndices,...
    shrinkParams,phaseErrorScale,shrinkErrorScale,verbose);
myObj = @(pRGC) myFit(pRGC,false);
myNonbcon = @(pRGC) nonbcon(pRGC,nBlockParams,nEccBands,nCellClasses);


%% Options
% The objective function is deterministic
optionsBADS.UncertaintyHandling = 0;


%% Search
if rgcSearchFlag
    [pRGC, fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);
else
    pRGC = p0;
    fVal = myObj(pRGC);
end


%% Report results
% Call the objective at the solution to report the fVals
myFit(pRGC,verboseFlag);

% Print the parameters in a format to be used as a seed in future searches
if verboseFlag
    fprintf('%% fVal = %2.1f \n',fVal);
    str = 'p0 = [... \n ';
    for ss=1:length(pRGC)
        str = [str sprintf('%2.10f, ',pRGC(ss))];
        if mod(ss,nBlockParams)==0
            str = [str '...\n'];
        end
    end
    str = [str(1:end-2) ' ];\n'];
    fprintf(str);
end

% Reshape and store the parameters
LMRatio = pRGC(end);
coneDelay = pRGC(end-1);
cfCone = pRGC(end-2);
pRGC = reshape(pRGC(1:nBlockParams*nEccBands*nCellClasses),[nBlockParams,nEccBands,nCellClasses]);

% Report the common params
if verboseFlag
    fprintf('cfCone: %2.2f, coneDelay: %2.2f, LMRatio: %2.2f \n',cfCone,coneDelay,LMRatio)
end

%% Interpolate params across eccentricity

% This is the position within each RGC recording bin that the model has
% identified as the correct location
eccDegs = zeros(1,nEccBands);
for ee = 1:nEccBands
    eccDegs(ee) = eccBins{ee}(1)+squeeze(pRGC(nBlockParams,ee,1)*range(eccBins{ee}));
end

% A simple linear interpolation and extrapolation, bounded by maximum and
% minimum returned values from the search
myInterpObj = @(v,xq,ii) max([repmat(min(v),1,length(xq)); min([repmat(max(v),1,length(xq)); interp1(eccDegs,v,xq,'linear','extrap')])]);

% Loop across cells
for cc=1:nCellClasses
    % Loop across the 7 params that vary with eccentricty
    for ii=1:nBlockParams-1
        % The values for this param across eccentricity, and the mean
        y = squeeze(pRGC(ii,:,cc));
        pMean(ii,cc) = mean(y);
        % Fit the values and store the fit
        pFitByEccen{ii,cc} = @(xq) myInterpObj(y,xq,ii);
    end
end


%% Save the temporalModel
rgcTemporalModel.p = pRGC;
rgcTemporalModel.LMRatio = LMRatio;
rgcTemporalModel.coneDelay = coneDelay;
rgcTemporalModel.cfCone = cfCone;
rgcTemporalModel.pMean = pMean;
rgcTemporalModel.pFitByEccen = pFitByEccen;
rgcTemporalModel.meta.blockParamNames = blockParamNames;
rgcTemporalModel.meta.ubBlock = ubBlock;
rgcTemporalModel.meta.cellClassIndices = cellClassIndices;
rgcTemporalModel.meta.eccFields = eccFields;
rgcTemporalModel.meta.eccBins = eccBins;
rgcTemporalModel.meta.eccDegs = eccDegs;
rgcTemporalModel.meta.lbBlock = lbBlock;
rgcTemporalModel.meta.ubBlock = ubBlock;

savePath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
save(savePath,'rgcTemporalModel');

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions


%% nonbcon
function c = nonbcon(p,nBlockParams,nEccBands,nCellClasses)

% Sometimes BADS sends an empty set of parameters
if isempty(p)
    c=1;
end

% Enforce that some parameters, such as delay and filter frequency,
% increase in value across eccentricity
for ii=1:size(p,1)
    subP = p(ii,:);
    subP = reshape(subP(1:nBlockParams*nEccBands*nCellClasses),[nBlockParams,nEccBands,nCellClasses]);
    for cc=1:nCellClasses
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
if nCellClasses>1
    c=any(c')';
end

end

%                 any(diff(tempP(4,:))<0) || ... % force cf2ndStage to increase with eccentricity
%                 any(diff(tempP(5,:))<0) || ... % force 2nd stage Q to increase with eccentricity




%% rgcFitOverallObjective
function fVal = rgcFitOverallObjective(pRGC,rgcData,eccFields,eccBins,nBlockParams,cellClassIndices,shrinkParams,phaseErrorScale,shrinkErrorScale,verbose)
% Given the midget and parasol data, as well as the full set of parameters
% across cell classes and eccentricity, derive the model fit error, with
% separate initial terms for fitting the gain and phase components of the
% transfer function. Additionally, we implement a "shrink" error that is
% used to drive certain sets of parameters to be the same across
% eccentricity locations.

cellClassData = fields(rgcData);
nCellClassData = length(cellClassData);
nCellClassModel = length(cellClassIndices);
nEccBins = length(eccBins);

% Extract and reshape the parameters
LMRatio = pRGC(end);
coneDelay = pRGC(end-1);
cfCone = pRGC(end-2);
pRGC = reshape(pRGC(1:nBlockParams*nEccBins*nCellClassModel),[nBlockParams,nEccBins,nCellClassModel]);

% Loop through the cell classes
for cc = 1:nCellClassData

    % Which cell class are we working on
    thisCellClass = cellClassData{cc};

    % Figure out which block of parameters corresponds to this cell class
    cellClassIndex = strcmp(cellClassIndices,thisCellClass);

    % Get the parameter block for this class
    pRGCCellBlock = squeeze(pRGC(:,:,cellClassIndex));

    % Calculate the gain and phase errors in fitting these data
    [fValGain(cc),fValPhase(cc)]  = ...
        rgcFitObjectiveByClass(pRGCCellBlock,cfCone,coneDelay,LMRatio,rgcData,thisCellClass,eccFields,eccBins);

    % Apply the phaseErrorScale
    fValPhase(cc) = fValPhase(cc) * phaseErrorScale;

    % Calculate the shrink error for this class
    pSub = squeeze(pRGC(shrinkParams{cellClassIndex},:,cellClassIndex));
    fValShrink(cc) = shrinkErrorScale * norm(std(pSub,[],2)./mean(pSub,2));

    % Report the values
    if verbose
        fprintf([thisCellClass ' - fValGain: %2.2f, fValPhase: %2.2f, fValShrink: %2.2f \n'],fValGain(cc),fValPhase(cc),fValShrink(cc));
    end

end

% The overall error is the sum of the individual errors
fVal = sum(fValGain) + sum(fValPhase) + sum(fValShrink);

end


%% rgcFitObjectiveByClass
function [fValGain,fValPhase]  = rgcFitObjectiveByClass(pRGCCellBlock,cfCone,coneDelay,LMRatio,rgcData,cellClass,eccFields,eccBins)

% Loop across eccentricity bands.
for ee = 1:length(eccFields)

    % The eccentricity domain
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % Extract the parameter values for this eccentricity band
    pRGCBlock = squeeze(pRGCCellBlock(:,ee));

    % The eccProportion parameter is used to determine the precise
    % eccentricity location within the range provided by the eccBin.
    eccProportion = pRGCBlock(8);
    eccDeg = eccBin(1)+eccProportion*(range(eccBin));


    % Obtain the chromatic weights and set the stimulus directions.
    stimulusDirections = {};
    switch cellClass
        case 'midget'
            stimulusDirections = {'LminusM','LMS'};
        case 'parasol'
            stimulusDirections = {'LMS'};
        case 'bistratified'
            stimulusDirections = {'S'};
    end

    % Obtain the chromatic weights
    chromaticCenterWeight = []; chromaticSurroundWeight = [];
    for ii=1:length(stimulusDirections)
        [chromaticCenterWeight(ii),chromaticSurroundWeight(ii)] = ...
            returnRGCChromaticWeights(cellClass,stimulusDirections{ii},eccDeg,LMRatio);
    end

    % Loop through the stimulus directions
    for ss = 1:length(stimulusDirections)

        thisStimulusDirection = stimulusDirections{ss};

        % Get the data
        freqsData = rgcData.(cellClass).(eccField).(thisStimulusDirection).f;
        gainData = rgcData.(cellClass).(eccField).(thisStimulusDirection).g;
        if isfield(rgcData.(cellClass).(eccField).(thisStimulusDirection),'p')
            phaseData = rgcData.(cellClass).(eccField).(thisStimulusDirection).p;
        else
            phaseData = [];
        end

        % Get the temporal RF for this cell and direction, and evaluate
        % at the specified temporal frequencies
        rfRGC = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight(ss),chromaticSurroundWeight(ss));
        rgcTTF = double(subs(rfRGC,freqsData));

        % Obtain the error in fitting the gain
        fValGain(ee,ss) = norm(gainData - abs(rgcTTF));

        % If there are phaseData, obtain the error in fitting
        if ~isempty(phaseData)
            fValPhase(ee,ss) = norm(phaseData - unwrap(angle(rgcTTF))*(180/pi));
        else
            fValPhase(ee,ss) = 0;
        end

    end % Loop over stimulus directions

end % Loop over eccentricity

% Return the norm of the errors across eccentricity
fValGain = norm(fValGain(:));
fValPhase = norm(fValPhase(:));

end