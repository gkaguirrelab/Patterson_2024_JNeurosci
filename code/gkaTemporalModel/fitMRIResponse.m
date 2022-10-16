function [p,fVal] = fitMRIResponseData(p0, v1FreqX, v1Eccentricity, v1Y, v1W, lgnFreqX, lgnY, lgnW, whichModel, useMonotonicConstraint)
% Fit the RGC-referred temporal model to combined V1 and LGN data
%
% Syntax:
%   output = myFunc(input)
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
%   p0                    - 1x[4+k*2] vector of parameters for the model fit.
%   v1FreqX               - 1xn vector of frequency values (in Hz) the
%                           correspond to the cortical amplitude values in
%                           v1Y. The length is equal to j * k, where j is
%                           the number of frequencies and k is the number
%                           of eccentricities.
%   v1Eccentricity        - 1xn vector of eccentricity locations (in
%                           degrees) from which the v1Y measurements were
%                           made.
%   v1Y                   - 1xn vector of BOLD amplitudes at each of many
%                           eccentricities and stimulus frequencies. The
%                           order of these is:
%                             [ [e1f1 e1f2 ... e1fj] [e2f1 e2f2 ... e2fj]
%   v1W                   - 1xn vector of weights for each of the
%                           measurements in v1Y. One choice is the inverse
%                           of the 95% CI of the measurement.
%   lgnFreqX              - 1xj vector of frequencies for which the LGN
%                           measurements were made.
%   lgnY                  - 1xj vector of BOLD amplitudes measured from the
%                           LGN.
%   lgnW                  - 1xj vector of weights for the measurements in
%                           lgnY.
%   whichModel            - String or char vec from the set:
%                               {'chromatic','luminance'}
%   useMonotonicConstraint - Logical. Controls if the model includes a
%                           non-linear constraint that requires the
%                           surround index decrease in value across
%                           eccentricity positions.
%
% Outputs:
%   p                     - 1x16 vector of model fit parameters

% Extract this value for later
nEcc = length(unique(v1Eccentricity));

% Define the objective, plausible bound, and p0, depending upon subject and
% stimulus
switch whichModel

    case 'chromatic'

        % Returns the TTF, and handles reshaping into a linear vector
        myV1TTF = @(p) assembleV1ChromResponseAcrossEcc(p,v1FreqX,v1Eccentricity);
        myLGNTTF = @(p) returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity);

        % Plausible bounds for the search
        plb = [ 0.01 10 0.5 20 repmat(0.2,1,nEcc) zeros(1,nEcc)];
        pub = [ 0.10 20 1.0 30 repmat(0.8,1,nEcc) ones(1,nEcc)];

    case 'luminance'

        % Returns the TTF, and handles reshaping into a linear vector
        myV1TTF = @(p) assembleV1LumResponseAcrossEcc(p,v1FreqX,v1Eccentricity);
        myLGNTTF = @(p) returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity);

        % Plausible bounds for the search
        plb = [ 0.01 20 0.1 10 repmat(0.2,1,nEcc) zeros(1,nEcc)];
        pub = [ 0.10 30 0.3 20 repmat(0.8,1,nEcc) ones(1,nEcc)];

end

% The weighted objective
myObj = @(p) norm(v1W.*(v1Y - myV1TTF(p))) + ...
    norm(lgnW.*(lgnY - myLGNTTF(p)));

% hard bounds
lb = [ 0  10 0.01 05 zeros(1,nEcc) zeros(1,nEcc)];
ub = [ 1  50 2.00 40 ones(1,nEcc) ones(1,nEcc)];

% Non-linear constraint that surround index decreases with eccentricity
if useMonotonicConstraint
    myNonbcon = @(p) nonbcon(p);
else
    myNonbcon = [];
end

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;
optionsBADS.Display = 'iter';

% search
[p,fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

end % main function



%% LOCAL FUNCTIONS

%% nonbcon

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(p)
nEcc = 6; nFixed = 4;
surroundIndex = p(:,nFixed+1:nFixed+nEcc);
c = sum(diff(surroundIndex,1,2)>0,2);
end

function response = assembleV1ChromResponseAcrossEcc(p,v1FreqX,v1Eccentricity)
% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
% Info needed to unpack the param vector
nFixed = 4;
nEcc = length(eccDegVals);
% Build the response vector
response = [];
parfor ee=1:length(eccDegVals)
    pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
    response(ee,:) = returnV1ChromEccTTFFit(pBlock,studiedFreqs,eccDegVals(ee));
end
response = reshape(response',1,length(v1FreqX));
end


function response = assembleV1LumResponseAcrossEcc(p,v1FreqX,v1Eccentricity)
% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
% Info needed to unpack the param vector
nFixed = 4;
nEcc = length(eccDegVals);
% Build the response vector
response = [];
parfor ee=1:length(eccDegVals)
    pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
    response(ee,:) = returnV1LumEccTTFFit(pBlock,studiedFreqs,eccDegVals(ee));
end
response = reshape(response',1,length(v1FreqX));
end
