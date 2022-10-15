function p = fitMRIResponseData(p0, v1FreqX, v1Eccentricity, v1Y, v1W, lgnFreqX, lgnY, lgnW, whichModel, useMonotonicConstraint  )
% Fit the RGC-referred temporal model to combined V1 and LGN data
%
% Syntax:
%   output = myFunc(input)
%
% Description:
%   Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean euismod
%   nulla a tempor scelerisque. Maecenas et lobortis est. Donec et turpis
%   sem. Sed fringilla in metus ut malesuada. Pellentesque nec eros
%   efficitur, pellentesque nisl vel, dapibus felis. Morbi eu gravida enim.
%   Sed sodales ipsum eget finibus dapibus. Fusce sagittis felis id orci
%   egestas, non convallis neque porttitor. Proin ut mi augue. Cras posuere
%   diam at purus dignissim, vel vestibulum tellus ultrices
%
% Inputs:
%   none
%   foo                   - Scalar. Foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo
%
% Optional key/value pairs:
%   none
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
	foo = 1;
    bar = myFunc(foo);
	fprintf('Bar = %d \n',bar);   
%}

% Fit the Mt Sinai fMRI data with a model that starts with the RGC temporal
% sensitivity functions

% Plotting mode. Copy and paste this to the console
%{
    searchFlag = true;
    for whichStim = 1:3; for whichSub = 1:2; fitMtSinaiResponseData; end; end
%}

% Search mode. Copy and paste this to the console
%{
    searchFlag = true;
    whichStim = 1; whichSub = 1;
    fitMtSinaiResponseData
%}

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

% Non-linear constraint that surround index decrease with eccentricity
if useMonotonicConstraint
    myNonbcon = @(p) nonbcon(p);
else
    myNonbcon = [];
end

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;
optionsBADS.Display = 'iter';

% search
p = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

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
