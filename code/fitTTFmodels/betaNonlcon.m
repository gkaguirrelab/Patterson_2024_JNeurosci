function y = betaTTF(p,w)
% Fits a non-central beta to temporal sensitivity data
%
% Syntax:
%  y = betaTTF(p,w)
%
% Description:
%   The fitting function is a non-central beta, that is further modified to
%   allow adjustment of the bounded interval and scaling of the overall
%   amplitude. The function is constrained to hold the first degree
%   parameter to an arbitrarily small value, and to apply non-linear
%   constraints in the fitting (described below).
%
%   There is no particular theoretical motivation for using this fitting
%   form. The fit does reflect the following expectations:
%     - The amplitude of response will approach zero as the stimulus
%       frequency approaches 1 Hz.
%     - The amplitude of the response will return to zero at higher
%       frequencies
%
%   When fitting noisy data, it may be desirable to also enforce the
%   non-linear constraints provided by betaNonlcon
%
%   Modified from:
%       Joshua D Carmichael
%       josh.carmichael@gmail.com
%
%   This function computes the probability density function for the
%   noncentral beta distribution using a transformation of variables to put
%   the desired density function in terms of a noncentral F, which is
%   included in Matlabs statsitics toolbox already.
%
% Inputs:
%   p                     - 1x3 vector of parameter values:
%                            a: The overall amplitude of the response  
%                            b: scaling of frequency indicies
%                           s2: The second degree of freedom/shape param
%                           nc:	The noncentrality parameter.
%   w                     - 1xn vector of temporal frequencies in Hz
%
% Outputs:
%   y                     - 1xn vector of values that is the modeled
%                           response
%
% Examples:
%{
    % Basic functionality
    w = [2 4 8 16 32 64 128];
    p = [1 8 1 1];
    y = betaTTF(p,w);
    semilogx(w,y,'-r');
%}
%{
    % Use in a non-linear search. Here are the data and studied frequencies
    w = [2 4 8 16 32 64 128];
    y = [0.6918    1.5463    2.7830    4.8807    3.7902    0.1244 0];

    % Create a scaled-up, log-spaced, version of the frequency domain
    wDelta = min(diff(log10(w)));
    upScale = 10;
    wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

    % Set up the p0 guess, and the bounds on the params
    p0 = [1 8 1 1];
    lb = [0 0 0 1];
    ub = [100 8 100 100];

    % The options for the search (mostly silence diagnostics)
    options = optimoptions(@fmincon,...
        'Diagnostics','off',...
        'Display','off');

    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - betaTTF(p,w));

    % Search
    p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

    % Obtain the high-resolution model fit
    yFit = betaTTF(p,wFit);

    % Plot it
    figure
    semilogx(w,y,'xk');
    hold on
    semilogx(wFit,yFit,'-r');
%}

% Convert w (in freq) to wIdx (which is 2^x)
wIdx = log10(w)./log10(2);

% Decompose p to the individual params
a = p(1);
b = p(2);
s2 = p(3);
nc = p(4);

% Perform the calculation
const = 1e-6./s2;
yFunc = @(r) ncfpdf( r./(const*(1-r)), 1e-6, s2, nc).* 1./(const*(1-r).^2);
y = a.*yFunc(wIdx./b);



end

