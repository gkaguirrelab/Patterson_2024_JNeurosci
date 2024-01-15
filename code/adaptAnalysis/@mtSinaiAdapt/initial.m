function x0 = initial(obj)
% Returns initial guess for the model parameters
%
% Syntax:
%   x0 = obj.initial()
%
% Description:
%   Initial values for the prf_timeShift model. Rationale is as follows:
%       x, y :  Center of the stimulus
%       sigma:  1 or 10 pixels, depending upon obj.scale
%       gain :  Set by obj.typicalGain
%       exp  :  Locked to 0.05, following Benson et al, 2018, HCP 7T data
%       shift:  Zero HRF temporal shift      
%
% Inputs:
%   none
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   x0                    - 1xnParams vector.
%


% Obj variables
typicalGain = obj.typicalGain;
nParams = obj.nParams;
nGainParams = obj.nGainParams;
nAdaptParams = obj.nAdaptParams;

% Assign the x0 variable
x0 = zeros(1,nParams);

% x0 gain
x0(1:nGainParams) = typicalGain;

% x0 exponential adaptation (tau in seconds; asymptote in proportion max)
tauIdx = nGainParams+1:1:nGainParams+nAdaptParams;
x0(tauIdx) = 12;

% x0 HRF: Flobs population mean amplitudes
x0(nGainParams+nAdaptParams+1:nParams) = [0.86, 0.09, 0.01];


end

