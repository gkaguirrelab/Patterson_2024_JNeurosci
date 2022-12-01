function results = fitMRIResponse(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,modelType,useMonotonicConstraint,paramSearch,verbose)
% Fit the RGC-referred temporal model to combined V1 and LGN data
%
% Syntax:
%  results = fitMRIResponseData(p0,stimulusDirections,studiedEccentricites,studiedFreqs,v1Y,v1W,lgnY,lgnW,useMonotonicConstraint,modelType)
%
% Description:
%   Model that simultaneously fits TTF responses from the LGN, and from V1
%   across eccentricity for multiple post-receptoral directions. The model
%   starts with Fourier-domain models of the chromatic and achromatic
%   temporal sensitivity functions of the RGCs as a function of
%   eccentricity. The model imposes a low-pass filter upon the modeled
%   responses at the retino-geniculate synapse, and again at the
%   geniculo-cortico synapse. At the level of V1, responses for a
%   particular post-receptoral direction (achromatic, chromatic red-green,
%   chromatic blue-yellow) are pooled and subject to delayed surround
%   inhibition. The surround delay varies by post-receptoral direction, and
%   the surround weight varries by by both direction and eccentricity. The
%   fit of the responses to the fMRI data are controlled by gain parameters
%   at the LGN and V1 stages.
%
%   The effect of eccentricity is present in the model in both a fixed and
%   parameterized form. There is a fixed effect of eccentricity at the
%   retina, which is derived from the change in the center and surround
%   cone inputs to the midget RGCs, and by changes in the parameters of the
%   temporal model used to fit the empirical RGC data. Further, the model
%   accounts for the variation in RGC cell populations with retinal
%   eccentricity, and the change in annular retinal area with distance from
%   the fovea. The effect of eccentricity in the post-retinal model is
%   found in the varying of the index of surround inhibition, and the gain
%   parameter that scales the resulting temporal sensitivity function to
%   the observed BOLD fMRI data.
%
% Inputs:
%   p0                    - 1x[5+c+s+s*k*2] vector of parameters for the
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
%   modelType             - Char vector. Options are "cell" or "stimulus".
%                           Controls if the model parameters operate upon
%                           the separate cell classes, or upon the separate
%                           stimuli. This distinction matters for the gain
%                           parameters at the lgn and v1 level, and for the
%                           delayed surround subtraction at V1.
%   useMonotonicConstraint - Logical. Controls if the model includes a
%                           non-linear constraint that requires the
%                           surround index to decrease in value across
%                           eccentricity positions for a given cell class.
%   paramSearch           - Char vector. Defines which parameters will be
%                           free to vary in the search, and if the v1 or
%                           lgn objectives will be used. This is to allow
%                           model fitting in stages in the initial
%                           parameter search, as there are many parameters.
%
% Outputs:
%   results               - Structure, with the fields:
%                           - pMRI: 1x[2+c+6+s*k*2] vector fit parameters
%                           - fVal: the objective function at the solution
%                           - cellClasses: the input classes
%                           - stimulusDirections: the input directions
%                           - paramCounts: Information regarding the number
%                             of parameters of each type. Used to unpack
%                             the parameter vector.
%


% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Define a variable with the number of eccentricities
nEccs = length(studiedEccentricites);

% The three types of cell classes in the model
cellClasses = {'midget','bistratified','parasol'};

% Define the bounds
lb = []; plb = []; pub = []; ub = [];

% "Unique" params. These implement filters on the different cell classes,
% and a neural->BOLD non-linearity
% - Non-linearity of neural --> BOLD activity
% - retino-geniculate, second order filter corner freq
% - geniculo-striate, second order filter corner freq
lb =  [ lb 0.5 repmat(20,1,2)];
plb = [plb 0.8 repmat(21,1,2)];
pub = [pub 0.9 repmat(40,1,2)];
ub =  [ ub 1.0 repmat(50,1,2)];
paramCounts.unique = 3;

% LGN BOLD fMRI gain, organized by stimulus
% - surround index
% - BOLD response gain
% The surround delay is the same as that used at the V1 level.
for cc = 1:length(stimulusDirections)
    lb =  [ lb 0 0];
    plb = [plb 0.5 1];
    pub = [pub 0.8 10];
    ub =  [ ub 1 100];
end
paramCounts.lgn = 2;

% V1 params. These vary by stimulus direction
% - surround delay
% - surround index (varies with eccentricity)
% - BOLD response gain (varies with eccentricity)
for ss = 1:length(stimulusDirections)
    lb =  [ lb 3 repmat(-1,1,nEccs) zeros(1,nEccs)];
    plb = [plb 5 zeros(1,nEccs) ones(1,nEccs)];
    pub = [pub 25 repmat(0.8,1,nEccs) repmat(10,1,nEccs)];
    ub =  [ ub 30 ones(1,nEccs) repmat(100,1,nEccs)];
end

% Some info about the params we will use later to unpack the vector
paramCounts.v1fixed = 1;
paramCounts.v1eccen = nEccs*2;
paramCounts.v1total = paramCounts.v1fixed+paramCounts.v1eccen;

% The objectives, each of which returns TTFs, reshaped across stimulus
% directions and eccentricities into a single vector
myV1TTF = @(pMRI) assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,paramCounts,modelType);
myLGNTTF = @(pMRI) assembleLGNResponse(pMRI,cellClasses,stimulusDirections,studiedFreqs,rgcTemporalModel,paramCounts,modelType);

% Which parameters and objective(s) shall we use in the search? 
switch paramSearch
    case 'fullLGN'
        lockIdx = [1 5:7 17:54];
        myObj = @(pMRI) norm(lgnW.*(lgnY - myLGNTTF(pMRI)));
    case 'fullV1'
        lockIdx = [2:4, 8:16];
        myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI)));
    case 'full'
        lockIdx = [];
        myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
            norm(lgnW.*(lgnY - myLGNTTF(pMRI)));
end
lb(lockIdx) = p0(lockIdx); plb(lockIdx) = p0(lockIdx);
ub(lockIdx) = p0(lockIdx); pub(lockIdx) = p0(lockIdx);

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(pMRI) nonbcon(pMRI,cellClasses,paramCounts,nEccs);
else
    myNonbcon = [];
end

% Options. Indicate that the objective function is deterministic, and
% handle verbosity
optionsBADS.UncertaintyHandling = 0;
if verbose
    optionsBADS.Display = 'iter';
else
    optionsBADS.Display = 'off';
end

% search
[pMRI,fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

% assemble the results structure
results.pMRI = pMRI;
results.fVal = fVal;
results.paramCounts = paramCounts;
results.cellClasses = cellClasses;
results.modelType = modelType;
results.paramSearch = paramSearch;
results.stimulusDirections = stimulusDirections;
results.studiedEccentricites = studiedEccentricites;
results.studiedFreqs = studiedFreqs;
results.v1Y = v1Y;
results.v1W = v1W;
results.lgnY = lgnY;
results.lgnW = lgnW;

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
