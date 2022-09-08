function totalRGC = totalRGC( showPlots)
% Count functions for the sum of all retinal ganglion cells
%
% Syntax:
%  totalRGC = cell.totalRGC(showPlots)
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) as a function of
%   retinal eccentricity (in degrees).
%
% Examples:
%{
    totalRGC = cell.totalRGC(true)
%}

% Handle plotting
if nargin==0
    showPlots = false;
end


% The meridians over which the calculation is to be performed
cardinalMeridianAngles = [0 90 180 270];
cardinalMeridianNames = {'nasal' 'superior' 'temporal' 'inferior'};

%% Cell counts
% We call the rgcDisplacementMap functions which have a representation of
% the Curcio & Allen 1990 RGC density results. The rgcDisplacementMap
% toolbox is referenced in eccentricity units of degrees visual field,
% and provides densities in square degrees visual field.


% Loop over the specified meridians
for mm = 1:length(cardinalMeridianAngles)
    
    % Set up this meridian model element
    totalRGC(mm).label = cardinalMeridianNames(mm);
    totalRGC(mm).angle = cardinalMeridianAngles(mm);
    
    % Obtain a spline fit to the cell densities
    splineFit = getSplineFitToRGCDensitySqDegVisual(cardinalMeridianAngles(mm));

    % Nan optic disc points and save the anonymous function.
    totalRGC(mm).countsDegSq =  @(posDeg) ...        
        zeroOpticDiscPoints(splineFit(posDeg), posDeg, cardinalMeridianAngles(mm));

    if showPlots
        if mm == 1
            figure
        end
        plot(0:0.5:50,totalRGC(mm).countsDegSq(0:0.5:50));
        hold on
    end
    
end

end

