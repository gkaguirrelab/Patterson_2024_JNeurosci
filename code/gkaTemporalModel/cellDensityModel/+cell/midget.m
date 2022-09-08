function midget = midget( totalRGC, showPlots, midgetLinkingFuncParams )
% Size and count functions for the midget RGC class
%
% Syntax:
%  midget = cell.midget( totalRGC, showPlots )
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) as a function of
%   retinal eccentricity (in degrees).
%
% Examples:
%{
	totalRGC = cell.totalRGC();
	midget = cell.midget( totalRGC, true );
%}

if nargin == 1
    showPlots = false;
    midgetLinkingFuncParams = [6.8151   21.4511    0.4668    0.9577];
end

if nargin == 2
    midgetLinkingFuncParams = [6.8151   21.4511    0.4668    0.9577];
end



% Define a maximum eccentricity and support of the model
maxEccenDeg = 50;
supportDeg = 0:0.01:maxEccenDeg;


% Obtain the midget fraction logistic function
[~, logisticFunc] = daceyMidgetFractionByEccenDegVisual();
midgetFraction = logisticFunc(...
    midgetLinkingFuncParams(1),...
    midgetLinkingFuncParams(2),...
    midgetLinkingFuncParams(3),...
    midgetLinkingFuncParams(4),...
    supportDeg);


%% Cell counts
% The midgetLinkingFuncParams are used to construct a function that
% converts the values within totalRGC into the midget fraction

% Define a support vector in visual degrees

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    % Obtain the total RGC density for this support
    countsDegSqTotal = totalRGC(mm).countsDegSq(supportDeg);
    
    % Set the optic disc nans to zero
    opticDiscPoints = isnan(countsDegSqTotal);
    
    % Obtain the midget counts
    countsDegSq = countsDegSqTotal .* midgetFraction;
    
    % Obtain a fit to the cell densities
    smoothFit = fit(supportDeg(~opticDiscPoints)', countsDegSq(~opticDiscPoints)', 'cubicinterp');
    
    % Set up this meridian model element
    midget(mm).label = totalRGC(mm).label;
    midget(mm).angle = totalRGC(mm).angle;
    
    % Nan optic disc points and save the anonymous function
    midget(mm).countsDegSq = @(posDeg) ...
        zeroOpticDiscPoints(smoothFit(posDeg), posDeg, totalRGC(mm).angle)';
    
    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        subplot(ceil(length(totalRGC)/2),ceil(length(totalRGC)/2),mm)
        plot(supportDeg,countsDegSqTotal,'-k');
        hold on
        plot(supportDeg,midget(mm).countsDegSq(supportDeg),'-r');
        title(totalRGC(mm).label);
    end
    
end


end


