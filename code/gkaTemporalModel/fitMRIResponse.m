function results = fitMRIResponse(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,useMonotonicConstraint,modelType)
% Fit the RGC-referred temporal model to combined V1 and LGN data
%
% Syntax:
%  [p,fVal] = fitMRIResponseData(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,useMonotonicConstraint)
%
% Description:
%   Model that simultaneously fits TTF responses from the LGN, and from V1
%   across eccentricity for multiple post-receptoral directions. The model
%   starts with Fourier-domain models of the chromatic and achromatic
%   temporal sensitivity functions of the RGCs as a function of
%   eccentricity. Responses at the LGN are modeled as the RGC responses,
%   subject to a delayed surround inhibition. Responses at V1 are the
%   responses from the LGN, subject to another iteration of delayed
%   surround inhibition, and also subject to a second-order, low-pass
%   temporal filter. Different sets of parameters modify the chromatic and
%   achromatic signals arriving from the LGN.
%
%   The effect of eccentricity is present in the model in both a fixed and
%   parameterized form. There is a fixed effect of eccentricity at the
%   retina, which is derived from the change in the center and surround
%   cone inputs to the midget RGCs, and by changes in the parameters of the
%   temporal model used to fit the empirical RGC data. The effect of
%   eccentricity in the post-retinal model is found in the varying of the
%   index of surround inhibition, and the gain parameter that scales the
%   resulting temporal sensitivity function to the observed BOLD fMRI data.
%
% Inputs:
%   p0                    - 1x[2+c+6+s*k*2] vector of parameters for the
%                           model fit, where c is the number of cell
%                           classes (midget, parasol, bistratified), s is
%                           the number of stimulus directions (LminusM, S,
%                           LMS), and k is the number of eccentricities.
%   stimulusDirections    - 1xs cell array naming the stimulus directions,
%                           such as {'LminusM','S','LMS'}
%   studiedEccentricites  - 1xk array of the centers of the studied
%                           eccentricity bins.
%   studiedFreqs          - 1xn vector of frequency values (in Hz).
%   v1Y                   - 1x(s*k*n) vector of BOLD amplitudes at each of
%                           many eccentricities and stimulus frequencies.
%                           The order of these is:
%                             [ [s1e1f1 s1e1f2 ... s1e1fj] [s1e2f1 s1e2f2 ... s1e2fn]
%   v1W                   - 1x(s*k*n) vector of weights for each of the
%                           measurements in v1Y. One choice is the inverse
%                           of the 95% CI of the measurement.
%   lgnY                  - 1x(s*n) vector of BOLD amplitudes measured from
%                           the LGN.
%   lgnW                  - 1x(s*n) vector of weights for the measurements
%                           in lgnY.
%   useMonotonicConstraint - Logical. Controls if the model includes a
%                           non-linear constraint that requires the
%                           surround index to decrease in value across
%                           eccentricity positions for a given cell class.
%
% Outputs:
%   p                     - 1x[2+c+6+s*k*2] vector of model fit parameters

if nargin<10
    modelType = 'full';
end

shrinkScaleFactor = 10;

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract this value for later
nEccs = length(studiedEccentricites);

% The three types of cell classes in the model
cellClasses = {'midget','bistratified','parasol'};

% Define the bounds
lb = []; plb = []; pub = []; ub = [];

% "Unique" params. These do not vary by cell class, channel, or
% eccentricity:
% - Relative stimulus power of LMS and L-M  on midget cells
% - lgn second order filter corner freq
% - lgn second order flter "quality"
% - V1 second order filter corner freq
% - V1 second order flter "quality"
lb =  [ lb 0.9 70 0.5 10 0.1];
plb = [plb 0.9 75 0.5 15 0.4];
pub = [pub 1.0 80 0.6 20 0.5];
ub =  [ ub 1.1 85 0.7 30 0.7];
paramCounts.unique = 5;

% LGN params. These vary by cell class
% - BOLD response gain
for cc = 1:length(cellClasses)
    lb =  [ lb 0];
    plb = [plb 1];
    pub = [pub 10];
    ub =  [ ub 100];
end
paramCounts.lgn = 1;

% V1 params. These vary by stimulus direction
% - surround delay
% - surround index (varies with eccentricity)
% - BOLD response gain
for ss = 1:length(stimulusDirections)
    lb =  [ lb 5 zeros(1,nEccs) zeros(1,nEccs)];
    plb = [plb 10 repmat(0.2,1,nEccs) ones(1,nEccs)];
    pub = [pub 25 repmat(0.8,1,nEccs) repmat(10,1,nEccs)];
    ub =  [ ub 30 ones(1,nEccs) repmat(100,1,nEccs)];
end

% Some info about the params we will use later to unpack the vector
paramCounts.v1fixed = 1;
paramCounts.v1eccen = nEccs*2;
paramCounts.v1total = paramCounts.v1fixed+paramCounts.v1eccen;

% For  reduced models, we lock some parameters
switch modelType
    case 'redGreenOnly'
        lockIdx = [1 2 3 6:8 22:47];
        v1W(37:end)=0;
    case 'blueYellowOnly'
        lockIdx = [1:5 6:8 9:21 35:47];
        v1W(1:36,73:end)=0;
    case 'chromaticOnly'
        v1W(73:end)=0;
    case 'achromaticOnly'
        lockIdx = [1:5 6:8 9:34];
        v1W(1:72)=0;
    case 'lockMidget'
        lockIdx = [2 3 4:6 9 12 13:25];
        v1W(1:36)=0;
    case 'lockBistratified'
        lockIdx = [6 7:9 12 26:38];
        v1W(37:72)=0;
    case 'v1IndexGainOnly'
    case 'lgn'
        lockIdx = [1:5 9:47];
    case 'full'
        lockIdx = [1];
end
lb(lockIdx) = p0(lockIdx); plb(lockIdx) = p0(lockIdx);
ub(lockIdx) = p0(lockIdx); pub(lockIdx) = p0(lockIdx);

% Returns the TTF, and handles reshaping into a linear vector
myV1TTF = @(pMRI) assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,paramCounts);
myLGNTTF = @(pMRI) assembleLGNResponse(pMRI,cellClasses,stimulusDirections,studiedFreqs,rgcTemporalModel,paramCounts);

% Define a shrink penalty that attempts to equate params
myShrinkPenalty = @(pMRI) calculateShrinkPenalty(pMRI,shrinkScaleFactor,nFixedParams,nEccs);

% The weighted objective
myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
      norm(lgnW.*(lgnY - myLGNTTF(pMRI)));

%myObj = @(pMRI) norm(lgnW.*(lgnY - myLGNTTF(pMRI)));

%myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI)));

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(pMRI) nonbcon(pMRI,cellClasses,paramCounts,nEccs);
else
    myNonbcon = [];
end

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;
optionsBADS.Display = 'iter';

% search
[pMRI,fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

% assemble the results structure
results.pMRI = pMRI;
results.fVal = fVal;
results.cellClasses = cellClasses;
results.stimulusDirections = stimulusDirections;
results.paramCounts = paramCounts;

end % main function


%% LOCAL FUNCTIONS

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(pMRI,cellClasses,paramCounts,nEccs)

for ss=1:3
    indexStart = paramCounts.unique + paramCounts.lgn*length(cellClasses) + (ss-1)*paramCounts.v1total + paramCounts.v1fixed;
    paramIndices = indexStart+1:indexStart+nEccs;
    surroundIndex = pMRI(:,paramIndices);
    c(:,ss) = sum(diff(surroundIndex,1,2)>0,2);
end
c = sum(c,2);
end


% Shrink penalty encourages the LGN stage to have the same temporal delay
% and index of suppression across
function shrinkPenalty = calculateShrinkPenalty(pMRI,shrinkScaleFactor,paramCounts,nEccs)
nParamsPerCellBlock = nFixedParams+nEccs*3;
% for pp = 1:3
%     indexStart = (pp-1)*nParamsPerCellBlock+nFixedParams;
%     indexToShrinkSets{pp} = indexStart+1:indexStart+nEccs;
% end
% For these parameters, try to match them across parameter blocks
blockIndicesToMatch = [12:18];
for ii=1:length(blockIndicesToMatch)
    indexToShrinkSets{ii}=blockIndicesToMatch(ii):nParamsPerCellBlock:nParamsPerCellBlock*2+blockIndicesToMatch(ii);
end
shrinkPenalty = 0;
for ii=1:length(indexToShrinkSets)
    shrinkPenalty = shrinkPenalty + std(pMRI(indexToShrinkSets{ii}))./mean(pMRI(indexToShrinkSets{ii}));
end
shrinkPenalty = shrinkPenalty * shrinkScaleFactor;
end