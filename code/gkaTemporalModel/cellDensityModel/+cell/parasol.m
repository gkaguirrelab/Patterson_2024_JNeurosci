function parasol = parasol( totalRGC, midget, bistratified, showPlots )
% Size and count functions for the parasol RGC class
%
% Syntax:
%  parasol = cell.parasol( totalRGC, midget, bistratified )
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) and cell
%   diameter (in mm) as a function of retinal eccentricity (in degrees).
%
% Examples:
%{
    totalRGC = cell.totalRGC_Curcio();
    midget = cell.midget_fixed( totalRGC );
    bistratified = cell.bistratified(totalRGC);
    parasol = cell.parasol( totalRGC, midget, bistratified );
%}

if nargin == 3
    showPlots = false;
end


% Define a maximum eccentricity of the model
maxEccenDeg = 50;


%% Cell counts
% We have solid information of total RGC density and an estimate of midget
% RGC density. The only parasol density measurements are from macaque. We
% therefore infer the parasol density by assuming that:
%   densityParasol = densityTotalRGC - densityMidget - densityBistratified

% Define a support vector in visual degrees
supportDeg = 0:0.01:maxEccenDeg;

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    % Set up this meridian model element
    parasol(mm).label = totalRGC(mm).label;
    parasol(mm).angle = totalRGC(mm).angle;
    
    % Create a function for parasol density
    parasol(mm).countsDegSq = @(x) ...
        (totalRGC(mm).countsDegSq(x)- ...
        midget(mm).countsDegSq(x) - ...
        bistratified(mm).countsDegSq(x));
    
    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        subplot(ceil(length(totalRGC)/2),ceil(length(totalRGC)/2),mm)
        plot(supportDeg,totalRGC(mm).countsDegSq(supportDeg),'--k');
        hold on
        plot(supportDeg,parasol(mm).countsDegSq(supportDeg),'-k');
        plot(supportDeg,midget(mm).countsDegSq(supportDeg),'-r');
        title(totalRGC(mm).label);
        ylim([0 2500]);
        drawnow
    end
    
end


end