function [returnVal, ceq] = betaTTF(p,w,y)
% Fits a non-central beta to temporal sensitivity data
%
% Syntax:
%  returnVal = betaTTF(p,w)
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
%   non-linear constraints. This provides for smoothness and prevent the
%   fit from greatly exceeding the range of the data.
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
%   y                     - 1xn vector data. Used to compute the non-linear
%                           constraint (optional). If passed, the routine
%                           returns the c constraint for returnVal.
%
% Outputs:
%   returnVal             - Either a 1xn vector of modeled values, or the
%                           scalar c value of the non-linear constraint
%   ceq                   - scalar ceq non-linear constraint value
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

    % The non-linear constraint
    myNonlcon = @(p) betaTTF(p,w,y);

    % Search
    p = fmincon(myObj,p0,[],[],[],[],lb,ub,myNonlcon,options);

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
yFit = a.*yFunc(wIdx./b);

% If no data were provided return the modeled values
if nargin==2
    returnVal = yFit;
    ceq = [];
else
    
    % Calculate the non-linear constraints. The rate of change should not
    % exceed 1 unit
    d = diff(yFit)./max(yFit);
    d = d(isfinite(d));
    c = max(abs(d))-1;
    
    % The interpolated peak should not be more than 10% the size of the
    % largest y value
    peak = (max(yFit)-max(y))./max(y) - 0.1;
    if peak > 0
        ceq = peak;
    else
        ceq = 0;
    end
    
    % Assign c to returnVal
    returnVal = c;
    
end

end

