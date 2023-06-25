function [p,fVal,yFit,yFitInterp] = fitWatsonModel(Y,W,studiedFreqs,interpFreqs)

% Handle inputs
if nargin == 3
    interpFreqs = logspace(log10(1),log10(100),501);
end

% fmincon options
options = optimoptions('fmincon');
options.Display = 'none';

% Set some bounds
LB = [0 1 0.5 0.5];
UB = [5 10 3 3];

% Two different p0 options
p0A = [1.5 5 1.1 1.5];
p0B = [4 1.5 1.5 1];

% Set up the objective
myObj = @(p) conditionedObj(p,Y,W,studiedFreqs,interpFreqs);

% Fit with two different p0 values
[pA, fValA] = fmincon(myObj,p0A,[],[],[],[],LB,UB,[],options);
[pB, fValB] = fmincon(myObj,p0B,[],[],[],[],LB,UB,[],options);
if fValA < fValB
    p = pA;
    fVal = fValA;
else
    p = pB;
    fVal = fValB;
end

% Get the values at the solution
yFit = watsonTemporalModel(p,studiedFreqs);
yFitInterp = watsonTemporalModel(p,interpFreqs);

end


function fVal = conditionedObj(p,Y,W,studiedFreqs,interpFreqs)

% Get the L2 norm of the weighted fit to the data
yFit = watsonTemporalModel(p,studiedFreqs);
fVal = norm( W .* ( Y - yFit) );

% Get the interpolated fit
yFitInterp = watsonTemporalModel(p,interpFreqs);

% Penalize responses that are not unimodal; this includes penalizing
% responses that have their highest value at either end of the frequency
% range
fVal = fVal + sum(sum(sign(diff(sign(diff(yFitInterp))))) == 0)*1e3;

end



