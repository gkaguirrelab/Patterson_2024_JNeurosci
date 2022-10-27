function [pMRI,fVal] = fitMRIResponse(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,useMonotonicConstraint)
% Fit the RGC-referred temporal model to combined V1 and LGN data
%
% Syntax:
%  [p,fVal] = fitMRIResponseData(p0, v1FreqX, v1Eccentricity, v1Y, v1W, lgnFreqX, lgnY, lgnW, whichModel, useMonotonicConstraint)
%
% Description:
%   Model that simultaneously fits TTF responses from the LGN, and from V1
%   across eccentricity for a particular post-receptoral direction. The
%   model starts with Fourier-domain models of the chromatic and achromatic
%   temporal sensitivity functions of the RGCs as a function of
%   eccentricity. Responses at the LGN are modeled as the RGC responses,
%   subject to a delayed surround inhibition. Responses at V1 are the
%   responses from the LGN, subject to another iteration of delayed
%   surround inhibition, and also subject to a second-order, low-pass
%   temporal filter.
%
%   The effect of eccentricity is present in the model in
%   both a fixed and parameterized form. There is a fixed effect of
%   eccentricity at the retina, which is derived from the change in the
%   center and surround cone inputs to the midget RGCs, and by changes in
%   the parameters of the temporal model used to fit the empirical RGC
%   data. The effect of eccentricity in the post-retinal model is found in
%   the varying of the index of surround inhibition.
%
%   There is one unique parameter of the model, shared across all cell
%   types and eccentricites:
%       LMConeRatio   - Ratio of L-to-M cones in the modeled retina
%
%   There are four parameters that are "fixed" for each cell type across
%   eccentricities:
%
%       lgnGain       - Multiplicative gain for fitting the LGN amplitudes
%       secondOrderFc - Corner frequency of 2nd order filter at V1 stage
%       secondOrderQ  - "Quality" parameter of 2nd order filter at V1
%       surroundDelay - The delay (in msecs) between the center and
%                       surround, used at both the LGN and V1 stages
%
%   and, there are two parameters that "float" across both cell type and
%   eccentricity:
%       surroundIndex - The index (range 0 - 1) of the surround inhibition
%                       at the LGN and V1 stages. One value for each of the
%                       k eccentricity measurements from V1
%       v1gain        - Multiplicative gain for fitting the V1 amplitudes
%                       at each eccentricity band.
%
% Inputs:
%   p0                    - 1x[c*(4+k*2)] vector of parameters for the
%                           model fit, where c is the number of cell
%                           classes (midget, parasol, bistratified)
%   v1FreqX               - 1x(s*n) vector of frequency values (in Hz) the
%                           correspond to the cortical amplitude values in
%                           v1Y. The length is equal to s*j*k, where s is
%                           the number of stimulusDirections
%                           ('LminusM','S','LMS'), j is the number of
%                           frequencies and k is the number of
%                           eccentricities.
%   v1Eccentricity        - 1x(s*n) vector of eccentricity locations (in
%                           degrees) from which the v1Y measurements were
%                           made.
%   v1Y                   - 1x(s*n) vector of BOLD amplitudes at each of
%                           many eccentricities and stimulus frequencies.
%                           The order of these is:
%                             [ [s1e1f1 s1e1f2 ... s1e1fj] [s1e2f1 s1e2f2 ... s1e2fj]
%   v1W                   - 1x(s*n) vector of weights for each of the
%                           measurements in v1Y. One choice is the inverse
%                           of the 95% CI of the measurement.
%   lgnFreqX              - 1x(s*j) vector of frequencies for which the LGN
%                           measurements were made.
%   lgnY                  - 1x(s*j) vector of BOLD amplitudes measured from
%                           the LGN.
%   lgnW                  - 1x(s*j) vector of weights for the measurements
%                           in lgnY.
%   stimulusDirections    - Cell array naming the stimulus directions, such
%                           as {'LminusM','S','LMS'}
%   useMonotonicConstraint - Logical. Controls if the model includes a
%                           non-linear constraint that requires the
%                           surround index to decrease in value across
%                           eccentricity positions for a given cell class.
%
% Outputs:
%   p                     - 1x[c*(4+k*2)] vector of model fit parameters

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract this value for later
nEcc = length(studiedEccentricites);
nUniqueParams = 9;
nFixedParams = 3;

% The influence of the shrink penalty
shrinkScaleFactor = 5;

% The model includes parameters for each of the cell classes
postReceptoralPaths = {'midget.LminusM','parasol.LMS','bistratified.S','midget.LMS'};

% Set up the parameters. The first set of "lgn" parameters are:
% - lgn surround delay (shared by all cell classes and eccentricity)
% - lgn surround index (shared by all cell classes and eccentricity)
% - lgn surround gain (varies by cell class)
% We then have blocks of parameters for each post-receptoral channel. First
% are the "fixed" parameters that do not vary by eccentricity:
% - second order filter corner frequency (Hz)
% - second order filter "quality" index (a.u.)
% - surround delay (msecs)
% The final set of parameters vary by post-receptoral pathway and by
% V1 eccentricity band:
% - surround index (a.u.)
% - surround gain (a.u.)

% lgn parameters, organized as midget, bistratified, parasol
lb =  repmat([5 0.1 0000],1,3);
plb = repmat([10 0.3 0.01],1,3);
pub = repmat([15 0.6 0.10],1,3);
ub =  repmat([20 0.7 1.00],1,3);

% v1 parameters, organized as LminusM, S, LMS-parasol, LMS-midget
for cc = 1:length(postReceptoralPaths)
    lb =  [ lb 10 0.1 10 zeros(1,nEcc) repmat(0.3,1,nEcc)];
    plb = [plb 15 0.3 12 repmat(0.2,1,nEcc) repmat(0.5,1,nEcc)];
    pub = [pub 25 0.6 30 repmat(0.8,1,nEcc) repmat(5,1,nEcc)];
    ub =  [ ub 90 0.7 40 ones(1,nEcc) repmat(100,1,nEcc)];
end

% Returns the TTF, and handles reshaping into a linear vector
myV1TTF = @(pMRI) assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nUniqueParams,nFixedParams);
myLGNTTF = @(pMRI) assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedFreqs,rgcTemporalModel);

% Shrink penalty that attempts to match various parameters
myShrinkPenalty = @(pMRI) calculateShrinkPenalty(pMRI,shrinkScaleFactor);

% The weighted objective
myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
    norm(lgnW.*(lgnY - myLGNTTF(pMRI))) + ...
    myShrinkPenalty(pMRI);

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(pMRI) nonbcon(pMRI,studiedEccentricites,postReceptoralPaths,nUniqueParams,nFixedParams);
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
function c = nonbcon(pMRI,studiedEccentricites,postReceptoralPaths,nUniqueParams,nFixedParams)
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
function shrinkPenalty = calculateShrinkPenalty(pMRI,shrinkScaleFactor)

indexToShrinkSets = {...
    [1 4 7] ... % The surround delay at the LGN stage
    [2 5 8] ... % The surround index at the LGN stage
    [12 27] ... % The V1 chromatic surround delays
    [42 57] ... % The V1 achromatic surround delays
    [10 25 40] ... % 2nd order V1 corner frequency for all except achromatic midget
%    [11 26 41] ... % 2nd order quality index for all except achromatic midget
    };

shrinkPenalty = 0;
for ii=1:length(indexToShrinkSets)
    shrinkPenalty = shrinkPenalty + std(pMRI(indexToShrinkSets{ii}))./mean(pMRI(indexToShrinkSets{ii}));
end

shrinkPenalty = shrinkPenalty * shrinkScaleFactor;

end
