
%% Housekeeping
clear
close all
verbose = true;


%% Pick a subject, stimulus, and eccentricity band
whichSub = 1;
whichStim = 1;
whichEcc = 1;


%% Load the data and basic model elements

% Load the Mt. Sinai data
mriData = loadMRIResponseData();

% Create the RGC temporal sensitivity model
rgcTemporalModel = fitRGCFResponse(false,false);

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
plotColor = {'r','b','k'};

% The number of acquisitions obtained for each measurement
nAcqs = 12;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Get the data
thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(['v1_ecc' num2str(whichEcc)]);
Y = mean(thisMatrix);
W = 1./std(thisMatrix);

% Define the response function
myResponse = @(p) returnResponse(p,stimulusDirections{whichStim},studiedEccentricites(whichEcc),studiedFreqs,rgcTemporalModel);

% Define the objective
myObj = @(p) norm(W.*(Y - myResponse(p)));

% No non-linear constraint at the moment
myNonbcon = [];


%% Set the bounds for the search
% The parameters are:
% - corner frequency [Hz] of synaptic filter retino-geniculate, and geniculo-v1
% - non-linear exponent for neural -> BOLD transform
% - gain of the response
% - delay (in msecs) for the surround
% - non-linear exponent for the divisive normalization
% - sigma for the divisive normalization
lb =  [  10 0.7   0   1 1.0 0.00];
plb = [  30 0.8   1  10 1.0 0.1];
pub = [  50 0.9  10  20 1.3 1] ;
ub =  [ 120 1.0 100 100 1.5 10];

% A p0 guess to start the search
p0 = [ 40    0.9    0.0687   30   1    1];

% Options. Indicate that the objective function is deterministic, and
% handle verbosity
optionsBADS.UncertaintyHandling = 0;
if verbose
    optionsBADS.Display = 'iter';
else
    optionsBADS.Display = 'off';
end

% The optimization toolbox is currently not available for Matlab running
% under Apple silicon. Detect this case and tell BADS so that it doesn't
% issue a warning
V = ver;
if ~any(strcmp({V.Name}, 'Optimization Toolbox'))
    optionsBADS.OptimToolbox = 0;
end

% BADS it
[p,fVal] = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);

% Plot it
figure
plot(Y,'-k');
hold on;

% Get the response
[yFit, rfFit] = myResponse(p);
plot(yFit,'*-r')

% Add a line for the best fit not using surround
pNoSurround = [119.9995    1.5000    0.0687         0    1    0];
plot(myResponse(pNoSurround),'-b')

% Bode plot of the temporal receptive field
plotRF(rfFit);



%% LOCAL FUNCTIONS

function [response, rfPostRetinal] = returnResponse(p,stimulusDirection,eccentricity,studiedFreqs,rgcTemporalModel)

secondOrderFc = p(1);
secondOrderQ = 0.45; % Fixed at 0.45
boldExp = p(2);
gain = p(3);
surroundDelay = p(4);
neuralExp = p(5);
sigma = p(6);

switch stimulusDirection
    case 'LminusM'
        activeCells = {'midget'};
    case 'S'
        activeCells = {'bistratified'};
    case 'LMS'
        activeCells = {'midget','parasol'};
end

for cc=1:length(activeCells)

    % Get the scaling effect of stimulus contrast
    stimulusContrastScale = returnStimulusContrastScale(activeCells{cc},stimulusDirection);

    % Get the post-retinal temporal RF
    rfPostRetinal(cc) = returnPostRetinalRF(...
        activeCells{cc},stimulusDirection,rgcTemporalModel,...
        eccentricity,stimulusContrastScale);

    % Second order low pass filter at the level of retino-
    % geniculate synapse
    rfPostRetinal(cc) = rfPostRetinal(cc).*stageSecondOrderLP(secondOrderFc,secondOrderQ);

    % Second order low pass filter at the level of genciulo-
    % cortical synapse
    rfPostRetinal(cc) = rfPostRetinal(cc).*stageSecondOrderLP(secondOrderFc,secondOrderQ);

    % Gain
    rfPostRetinal(cc) = (gain/1000)*rfPostRetinal(cc);

end

% Add the postRetinal RFs together
rfPostRetinal = sum(rfPostRetinal);

% If the delay variable is greater than zero, perform delayed surround
% divisive normalization
if surroundDelay > 0

    % An anonymous function to perform signed exponentiation
    signedExp = @(vec,k) sign(vec).*(abs(vec).^k);

    % The numerator for the divisive normalization
    numerator = signedExp(rfPostRetinal,neuralExp);

    % The denominator for the divisive normalization
    rfDelayed = rfPostRetinal.*stageDelay(surroundDelay/1000);
    denominator = signedExp(rfDelayed,neuralExp) + sigma^neuralExp;

    % Perform the division
    rfPostRetinal = numerator ./ denominator;

end

% Derive the amplitude and phase from the Fourier model
ttfComplex = double(subs(rfPostRetinal,studiedFreqs));
amplitude = abs(ttfComplex);
phase = unwrap(angle(ttfComplex));

% Apply the non-linear transformation of neural to BOLD response
response = amplitude.^boldExp;


end

