function [y,Hc,Hs] = watsonTTF2param(p,w)
% Implements the AB Watson center-surround temporal sensitivity model
%
% Syntax:
%  [y,Hc,Hs] = watsonTTF(p,w)
%
% Description
%   Watson citation and model details here
%
%   The model as implemented here fixes several parameters of the Watson
%   model based upon Figure 6.5 of Watson's Temporal Sensitivity chapter.
%   The degree of the center and surround filters are fixed at order 9 and
%   10, respectively, and the ratio of the surround-to-center time constant
%   is fixed at 1.33. This leaves the following free parameters:
%     - Overall amplitude (G)
%     - Time constant of the center filter (tc) 
%
% Inputs:
%   p                     - 1x2 vector of parameter values:
%   w                     - 1xn vector of temporal frequencies in Hz
%
% Outputs:
%   y                     - A 1xn vector of modeled values
%
% Examples:
%{
    % Basic functionality
    w = [2 4 8 16 32 64 128];
    p = [1.72, 0.87, 0.013];
    y = watsonTTF(p,w);
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
    p0 = [1.72, 0.87, 0.013];
    lb = [0 0 0];
    ub = [5 1 0.025];

    % The options for the search (mostly silence diagnostics)
    options = optimoptions(@fmincon,...
        'Diagnostics','off',...
        'Display','off');

    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - watsonTTF(p,w));

    % Search
    p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

    % Obtain the high-resolution model fit
    yFit = watsonTTF(p,wFit);

    % Plot it
    figure
    semilogx(w,y,'xk');
    hold on
    semilogx(wFit,yFit,'-r');
%}

% Decompose p to the individual params
G = p(1);
Gs = 0.9; % based on average Gs across subjects and directions V1 wide field
tc = p(2);
nc = 9;
ts = tc.*1.33;
ns = 10;

% Two linear bandpass filters
Hc=(w.*1i*2*pi*tc+1).^-1*nc;
Hs=(w.*1i*2*pi*ts+1).^-1*ns;

% Apply the magnitude scaling
Hs = Gs*Hs;

% The response is the difference between the center and surround
% filters, with the surround subject to a relative gain scaling
% (with a lower bound at zero) and the overall reesponse
y = G*(Hc-Hs);

% Return just the real component, as this is all we are currently using
y = real(y);

end
