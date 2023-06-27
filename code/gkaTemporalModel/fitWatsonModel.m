function [p,fVal,yFit,yFitInterp] = fitWatsonModel(Y,W,studiedFreqs,p0A,interpFreqs)

% Handle inputs
if nargin == 3
    p0A = [1.5 5 1.1 1.5];
    interpFreqs = logspace(log10(1),log10(100),501);
end

if nargin == 4
    interpFreqs = logspace(log10(1),log10(100),501);
end

% fmincon options
options = optimoptions('fmincon');
options.Display = 'none';

% Set some bounds
LB = [0 1 0.5 0.5];
UB = [10 10 3 3];

% Two different p0 options
p0B = [4 1.5 1.5 1];

% Set up the objective
myObj = @(p) objectiveFunc(p,Y,W,studiedFreqs);
myNonlcon = @(p) unimodalConstraint(p,interpFreqs);

% Fit with two different p0 values
[pA, fValA] = fmincon(myObj,p0A,[],[],[],[],LB,UB,myNonlcon,options);
[pB, fValB] = fmincon(myObj,p0B,[],[],[],[],LB,UB,myNonlcon,options);
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


function [c,ceq] = unimodalConstraint(p,interpFreqs)

% We don't need c
c = [];

% Get the interpolated fit
yFitInterp = watsonTemporalModel(p,interpFreqs);

% Penalize responses that are not unimodal
ceq = sum(sum(sign(diff(sign(diff(yFitInterp))))) == 0)*10;

end

function fVal = objectiveFunc(p,Y,W,studiedFreqs)

% Get the L2 norm of the weighted fit to the data
yFit = watsonTemporalModel(p,studiedFreqs);
fVal = norm( W .* ( Y - yFit) );

end



