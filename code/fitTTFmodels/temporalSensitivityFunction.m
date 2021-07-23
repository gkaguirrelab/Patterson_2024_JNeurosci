function y = temporalSensitivityFunction(p,w)
% Implements the AB Watson center-surround temporal sensitivity model
%
% Syntax:
%  y = temporalSensitivityFunction(p,w)
%
% Description
%   Watson citation and model details here
%
% Inputs:
%   p                     - 1x6 vector of parameter values:
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
    p = [1,0.3,0.02,1,1,3];
    y = temporalSensitivityFunction(p,w);
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
    p0 = [100,0.3,0.01,5,0.01,10];
    lb = [0 0 0 0 0 0];
    ub = [Inf 1 1 20 1 20];

    % The options for the search (mostly silence diagnostics)
    options = optimoptions(@fmincon,...
        'Diagnostics','off',...
        'Display','off');

    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - temporalSensitivityFunction(p,w));

    % Search
    p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

    % Obtain the high-resolution model fit
    yFit = temporalSensitivityFunction(p,wFit);

    % Plot it
    figure
    semilogx(w,y,'xk');
    hold on
    semilogx(wFit,yFit,'-r');
%}

% Decompose p to the individual params
G = p(1);
Gs = p(2);
tc = p(3);
nc = p(4);
ts = p(5);
ns = p(6);

% Two linear bandpass filters
Hc=(w.*1i*2*pi*tc+1).^-1*nc;
Hs=(w.*1i*2*pi*ts+1).^-1*ns;

% The response is the difference between the center and surround
% filters, with the surround subject to a relative gain scaling
% (with a lower bound at zero) and the overall reesponse
y=G*(Hc-(Gs*Hs));

% Return just the real component, as this is all we are currently using
y = real(y);

end
