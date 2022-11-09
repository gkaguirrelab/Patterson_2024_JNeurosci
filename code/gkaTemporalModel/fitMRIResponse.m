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
    modelType = 'full';
end

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract this value for later
nEcc = length(studiedEccentricites);
nUniqueParams = 9;
nFixedParams = 1;

% The model includes parameters for each of the "post receptoral paths",
% which is realized here as the interaction of retinal ganglion cell
% classes with stimulus directions.
postReceptoralPaths = {'midget.LminusM','bistratified.S','parasol.LMS','midget.LMS'};

% Set up the parameters.

% lgn parameters, organized as
% - lgn surround delay (msecs)
% - lgn surround index (proportion, in the range 0 - 1)
% - lgn gain (a.u.; varies by midget, bistratified, parasol)
lb =  [05 0.1 repmat(0000,1,3)];
plb = [10 0.3 repmat(0.01,1,3)];
pub = [15 0.6 repmat(0.10,1,3)];
ub =  [20 0.7 repmat(1.00,1,3)];

% v1 parameters organized as chromatic, achromatic
% - second order filter corner frequency (Hz)
% - second order filter "quality" index (a.u.)
for pp = 1:2
    lb =  [ lb 05 0.1];
    plb = [plb 15 0.3];
    pub = [pub 25 0.6];
    ub =  [ ub 90 0.7];
end

% v1 parameters that are organized by "pathway" that are fixed with
% eccentricity, organized as:
% LminusM, S, LMS-parasol, LMS-midget:
% - surround delay (msecs)
% v1 paramters that vary by eccentricity
% - surround index (a.u.)
% - surround gain (a.u.)
for pp = 1:length(postReceptoralPaths)
    lb =  [ lb 05 zeros(1,nEcc) repmat(0.3,1,nEcc)];
    plb = [plb 10 repmat(0.2,1,nEcc) repmat(0.5,1,nEcc)];
    pub = [pub 30 repmat(0.8,1,nEcc) repmat(5,1,nEcc)];
    ub =  [ ub 40 ones(1,nEcc) repmat(100,1,nEcc)];
end

% For  reduced models, we lock some parameters
switch modelType
    case 'rgc'
        % Zero out the surround weight
        lb([2 11:16 24:29 37:42 50:56]) = 0;
        % move the 2nd order filter out of range
        lb([6 8]) = 200; lb([7 9]) = 0.7;
        plb([2 11:16 24:29 37:42 50:56]) = 0;
        plb([6 8]) = 200; plb([7 9]) = 0.7;
        pub([2 11:16 24:29 37:42 50:56]) = 0;
        pub([6 8]) = 200; pub([7 9]) = 0.7;
        ub([2 11:16 24:29 37:42 50:56]) = 0;
        ub([6 8]) = 200; ub([7 9]) = 0.7;
        p0([2 11:16 24:29 37:42 50:56]) = 0;
        p0([6 8]) = 200; p0([7 9]) = 0.7;
    case 'lgn'
        lb([11:16 24:29 37:42 50:56]) = 0;
        lb([6 8]) = 200; lb([7 9]) = 0.7;
        plb([11:16 24:29 37:42 50:56]) = 0;
        plb([6 8]) = 200; plb([7 9]) = 0.7;
        pub([11:16 24:29 37:42 50:56]) = 0;
        pub([6 8]) = 200; pub([7 9]) = 0.7;
        ub([11:16 24:29 37:42 50:56]) = 0;
        ub([6 8]) = 200; ub([7 9]) = 0.7;
        p0([11:16 24:29 37:42 50:56]) = 0;
        p0([6 8]) = 200; p0([7 9]) = 0.7;
    case 'v1'
        lb([6 8]) = 200; lb([7 9]) = 0.7;
        plb([6 8]) = 200; plb([7 9]) = 0.7;
        pub([6 8]) = 200; pub([7 9]) = 0.7;
        ub([6 8]) = 200; ub([7 9]) = 0.7;
        p0([6 8]) = 200; p0([7 9]) = 0.7;
    case 'full'
        % Make no changes
end

% Returns the TTF, and handles reshaping into a linear vector
myV1TTF = @(pMRI) assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel,nUniqueParams,nFixedParams);
myLGNTTF = @(pMRI) assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedFreqs,rgcTemporalModel);
        
% The weighted objective
myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
    norm(lgnW.*(lgnY - myLGNTTF(pMRI)));

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
