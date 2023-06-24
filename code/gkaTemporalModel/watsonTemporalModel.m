function y = watsonTemporalModel(params, frequenciesToModel, filterOrder)
% Beau Watson's 1986 center-surround neural temporal sensitivity model
%
% Syntax:
%  y = watsonTemporalModel(params, frequenciesToModel, filterOrder)
%
% Description:
%	Calculates the two-component (center-surround) Watson temporal model
%   The parameters (p) defines both the "center" and the "surround" linear
%   filter components. The entire model is the difference between these two
%   filter components.
%
%	The model is expressed in Eq 45 of:
%
%       Watson, A.B. (1986). Temporal sensitivity. In Handbook of Perception
%       and Human Performance, Volume 1, K. Boff, L. Kaufman and
%       J. Thomas, eds. (New York: Wiley), pp. 6-1-6-43.
%
%   Note that there is a typo in original manuscript. Equation 45 should
%   be:
%
%       H(frequenciesHz) = a[H1(frequenciesHz) - bH2(frequenciesHz)]
%
%   where a and b are scale factors. We have modified the implementation of
%   scale factors here.
%
%   Additionally, Watson (1986) gives the time-constant of the model in
%   units of milliseconds, but we find that, to reproduce the presented
%   figures, this time-constant is converted at some point to units of
%   seconds prior to its entry into the equations.
%
%   The model is the difference between two linear impulse response
%   filters, each of which is themselves a cascade of low-pass filters. The
%   number of filters in the cascade (the "filterOrder") is set
%   empirically. Center and surround orders of "9" and "10" are presented
%   in (e.g.) Figure 6.5 of Watson (1986).
%
%   Note that the model can only return positive amplitudes of response.
%
% Inputs:
%   frequenciesToModel    - 1xn vector that provides the stimulus
%                           frequencies in Hz for which the model will be
%                           evaluated
%   params                - 1x4 vector of model parameters:
%                            gain - amplitude of the response at maximum
%                             tau - time constant of the center filter (in
%                                   msecs)
%                           kappa - multiplier of the surround time-constant
%                            zeta - relative amplitude of the surround
%   filterOrder           - 1x2 vector. Optional provides the order of the
%                           center and surround filters respetively. If not
%                           provided, defaults to [9 10]
%
% Outputs:
%   y                     - 1xn vector of modeled amplitude values.
%
% Examples:
%{
    % Demonstrate basic output of the model
    freqHz = logspace(0,log10(64),100);
    params = [3 2 2 1];
    y = watsonTemporalModel(params,freqHz);
    semilogx(freqHz,y,'-k');    
%}
%{
    % Fit the Watson model to some empirical data
    stimulusFreqHz = [0.5 1 2 4 8 16 32 64];
    pctBOLDresponse = [0    0.0670    0.2320    0.2850    0.5750    0.7160    0.4020    0.0030];

    % Prepare to search by defining an objective and setting bounds
    myObj = @(p) norm(pctBOLDresponse - watsonTemporalModel(p,stimulusFreqHz));
    LB = [0 1 0.5 0.5];
    UB = [5 10 3 3];
    p0 = [0.5 5 2 1.5];

    % Search
    p = fmincon(myObj,p0,[],[],[],[],LB,UB);

    % Plot the data and an interpolated fit
    stimulusFreqHzFine = logspace(-1,log10(100),100);
    semilogx(stimulusFreqHzFine,watsonTemporalModel(p,stimulusFreqHzFine),'-r');
    hold on
    semilogx(stimulusFreqHz, pctBOLDresponse, '*k');
    hold off
%}

if nargin == 3
    centerFilterOrder = filterOrder(1); % Order of the center (usually fast) filter
    surroundFilterOrder = filterOrder(2); % Order of the surround (usually slow) filter
else
    % Fixed parameters (taken from Figure 6.4 and 6.5 of Watson 1986)
    centerFilterOrder = 9; % Order of the center (usually fast) filter
    surroundFilterOrder = 10; % Order of the surround (usually slow) filter
end


% Define a frequency domain in Hz over which the model is defined. The
% maximum and minimum value of the y response should be contained within
% this range for any plausible set of parameters. We hard code a log-spaced
% range of 0.1 - 200 Hz here.
freqDomain = logspace(-1,log10(200),100);

% Sanity check the frequency input
if max(frequenciesToModel)>max(freqDomain) || min(frequenciesToModel)<min(freqDomain)
    error('The passed frequency is out of range for the model');
end

% Calculate the response. We need to search across values of the center
% amplitude to obtain a response that has a maximum value of unity
maxResponse = @(x) max(abs(returnModel(x, params, freqDomain, centerFilterOrder, surroundFilterOrder)));
myObj = @(x) abs(maxResponse(x) - 1);
params_centerAmplitude = fminsearch(myObj,1);

% The model calculates a vector of complex values that define the Fourier
% transform of the system output. As we are implementing a model of just
% the amplitude component of a temporal transfer function, the absolute
% value of the model output is returned.
y = params(1)*abs(returnModel(params_centerAmplitude, params, frequenciesToModel, centerFilterOrder, surroundFilterOrder));

end % main function


function y = returnModel(params_centerAmplitude, params, frequenciesToModel, centerFilterOrder, surroundFilterOrder)

% Un-pack the passed parameters
params_tau = params(2)./1000;   % Convert from msecs to secs
params_kappa = params(3);
params_zeta = params(4);

H1 = nStageLowPassFilter(params_tau,frequenciesToModel,centerFilterOrder);
H2 = nStageLowPassFilter(params_kappa*params_tau,frequenciesToModel,surroundFilterOrder);
y = (params_centerAmplitude * H1) - (params_zeta*params_centerAmplitude*H2);

end

function Hsub = nStageLowPassFilter(tau,frequenciesHz,filterOrder)
% This function implements the system response of the linear filter
% for temporal sensitivity of neural systems in Eq 42 of Watson (1986).
%
% The implemented function is the "system respone" (Fourier transform) of
% the impulse response of an nth-order filter which is of the form:
%
% h(t) = u(t) * (1/(tau*(n-1)!)) * (t/tau)^(n-1) * exp(-t/tau)
%
% tau -- Time constant of the filter (in seconds)
% frequenciesHz -- a vector of frequencies at which to realize the model
% filterOrder -- the number of low-pass filters which are cascaded in
%                the model
Hsub = (1i*2*pi*frequenciesHz*tau + 1) .^ (-filterOrder);

end