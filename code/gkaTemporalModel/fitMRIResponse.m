function [pMRI,fVal] = fitMRIResponse(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,useMonotonicConstraint,modelType)
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
    modelType = 'v1MidgetParasol';
end

shrinkScaleFactor = 1;

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract this value for later
nEccs = length(studiedEccentricites);

% The number of parameters that do not vary with eccentricity
nFixedParams = 5;

% The model includes parameters for each of the "post receptoral paths",
% which is realized here as the interaction of retinal ganglion cell
% classes with stimulus directions.
postReceptoralPaths = {'midget','bistratified','parasol'};

% Parameters are organized by cell "pathway". The first set are fixed with
% eccentricity;
% - lgn surround delay (msecs)
% - lgn surround index (proportion, in the range 0 - 1)
% - lgn gain (BOLD % / spikes / sec)
% - second order filter corner frequency (Hz)
% - second order filter "quality" index (a.u.)
%
% v1 paramters that vary by eccentricity
% - surround delay (msecs)
% - surround index (a.u.)
% - gain (a.u.)
lb = []; plb = []; pub = []; ub = [];
for pp = 1:length(postReceptoralPaths)
    % lgn parameters (surround delay, surround index, gain)
    lb =  [ lb 02 0.0 0000];
    plb = [plb 10 0.3 0.01];
    pub = [pub 15 0.7 0.10];
    ub =  [ ub 20 1.0 1.00];

    % V1 second order filter corner freq, quality
    lb =  [ lb 005 0.1];
    plb = [plb 015 0.3];
    pub = [pub 035 0.6];
    ub =  [ ub 100 0.7];

    % V1 params that vary by ecccentricity, surround delay, index, and
    % overall gain
    lb =  [ lb repmat(02,1,nEccs) zeros(1,nEccs) repmat(0,1,nEccs)];
    plb = [plb repmat(10,1,nEccs) repmat(0.2,1,nEccs) repmat(0.5,1,nEccs)];
    pub = [pub repmat(30,1,nEccs) repmat(0.8,1,nEccs) repmat(5,1,nEccs)];
    ub =  [ ub repmat(40,1,nEccs) ones(1,nEccs) repmat(100,1,nEccs)];
end

% A final parameter adjusts the relative stimulus effectiveness of LMS and
% L-M contrast upon midget cells.
lb  = [ lb 0.2];
plb = [plb 0.5];
pub = [pub 2.0];
ub  = [ ub 4.0];

% For  reduced models, we lock some parameters
nParamsPerCellBlock = nFixedParams+nEccs*3;
switch modelType
    case 'bistratifiedOnly'
        lockIdx = [1:nParamsPerCellBlock, nParamsPerCellBlock*2+1:nParamsPerCellBlock*3, 70];
    case 'v1MidgetOnly'
        lockIdx = [1+nParamsPerCellBlock:nParamsPerCellBlock*3, 70];
        v1W(37:end)=0;
    case 'v1ParasolOnly'
        lockIdx = [1:nParamsPerCellBlock*2, 70];
        v1W(1:end-36)=0;
    case 'v1MidgetParasol'
        lockIdx = [1+nParamsPerCellBlock:nParamsPerCellBlock*2];
        v1W(37:72)=0;
    case 'v1GainOnly'
        lockIdx = [1:17, 24:40, 47:63];
    case 'lgnOnly'
        lockIdx = [4:23, 27:46, 50:59];
    case 'full'
        lockIdx = [];
end
lb(lockIdx) = p0(lockIdx); plb(lockIdx) = p0(lockIdx);
ub(lockIdx) = p0(lockIdx); pub(lockIdx) = p0(lockIdx);

% Returns the TTF, and handles reshaping into a linear vector
myV1TTF = @(pMRI) assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nFixedParams);
myLGNTTF = @(pMRI) assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nFixedParams);
        
% Define a shrink penalty that attempts to equate params
myShrinkPenalty = @(pMRI) calculateShrinkPenalty(pMRI,shrinkScaleFactor,nFixedParams,nEccs);

% The weighted objective
% myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
%     norm(lgnW.*(lgnY - myLGNTTF(pMRI)));

myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
    myShrinkPenalty(pMRI);

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(pMRI) nonbcon(pMRI,studiedEccentricites,postReceptoralPaths,nFixedParams);
else
    myNonbcon = [];
end

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;
optionsBADS.Display = 'iter';

% search
[pMRI,fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

end % main function


%% LOCAL FUNCTIONS

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(pMRI,studiedEccentricites,postReceptoralPaths,nFixedParams)
nEccs = length(studiedEccentricites);
for whichCell=1:length(postReceptoralPaths)
    paramIndices = 1+nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams: ...
        nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams+nEccs;
    surroundIndex = pMRI(:,paramIndices);
    c(:,whichCell) = sum(diff(surroundIndex,1,2)>0,2);
end
c = sum(c,2);
end


% Shrink penalty encourages the LGN stage to have the same temporal delay
% and index of suppression across 
function shrinkPenalty = calculateShrinkPenalty(pMRI,shrinkScaleFactor,nFixedParams,nEccs)

nParamsPerCellBlock = nFixedParams+nEccs*3;

for pp = 1:3
    indexStart = (pp-1)*nParamsPerCellBlock+nFixedParams;
    indexToShrinkSets{pp} = [indexStart+1:indexStart+nEccs];
end

shrinkPenalty = 0;
for ii=1:length(indexToShrinkSets)
    shrinkPenalty = shrinkPenalty + std(pMRI(indexToShrinkSets{ii}))./mean(pMRI(indexToShrinkSets{ii}));
end

shrinkPenalty = shrinkPenalty * shrinkScaleFactor;

end