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
%   The parameters of the model are:
%       LMConeRatio   - Ratio of L-to-M cones in the modeled retina
%       lgnGain       - Multiplicative gain for fitting the LGN amplitudes
%       secondOrderFc - Corner frequency of 2nd order filter at V1 stage
%       secondOrderQ  - "Quality" parameter of 2nd order filter at V1
%       surroundDelay - The delay (in msecs) between the center and
%                       surround, used at both the LGN and V1 stages
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
nUniqueParams = 1;
nFixedParams = 4;

% The number of params that are fixed across eccentricity within each
% stimulus block


% The model includes parameters for each of the cell classes
cellClassOrder = {'midget','parasol','bistratified'};

% Set up the bounds. Initialize these vectors with the unique parameter
% LMRatio, that is commmon to all cell classes and eccentricities
lb = [0.33]; ub = [3]; plb = [0.5]; pub = [2];
for cc = 1:length(cellClassOrder)
    % hard bounds
    lb = [lb [ 0 10 0.01 05 zeros(1,nEcc) zeros(1,nEcc)]];
    ub = [ub [ 1 50 2.00 40 ones(1,nEcc) ones(1,nEcc)]];
    % Plausible bounds vary by stimulus direction
    switch cellClassOrder{cc}
        case 'midget'
            plb = [plb [ 0.01 10 0.5 20 repmat(0.2,1,nEcc) zeros(1,nEcc)]];
            pub = [pub [ 0.10 20 1.0 30 repmat(0.8,1,nEcc) ones(1,nEcc)]];
        case 'bistratified'
            plb = [plb [ 0.01 20 0.1 10 repmat(0.2,1,nEcc) zeros(1,nEcc)]];
            pub = [pub [ 0.10 30 0.3 20 repmat(0.8,1,nEcc) ones(1,nEcc)]];
        case 'parasol'
            plb = [plb [ 0.01 20 0.1 10 repmat(0.2,1,nEcc) zeros(1,nEcc)]];
            pub = [pub [ 0.10 30 0.3 20 repmat(0.8,1,nEcc) ones(1,nEcc)]];
    end
end

% Returns the TTF, and handles reshaping into a linear vector
myV1TTF = @(pMRI) assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams);
myLGNTTF = @(pMRI) assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,studiedFreqs,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams);

% The weighted objective
myObj = @(pMRI) norm(v1W.*(v1Y - myV1TTF(pMRI))) + ...
    norm(lgnW.*(lgnY - myLGNTTF(pMRI)));

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(pMRI) nonbcon(pMRI);
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

%% nonbcon

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(p)
nEcc = 6; nFixed = 5;
surroundIndex = p(:,nFixed+1:nFixed+nEcc);
c = sum(diff(surroundIndex,1,2)>0,2);
end