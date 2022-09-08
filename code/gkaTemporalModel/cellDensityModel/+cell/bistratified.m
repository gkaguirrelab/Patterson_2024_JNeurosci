function bistratified = bistratified( totalRGC, showPlots )
% Size and count functions for the bistratified RGC class
%
% Syntax:
%  bistratified = cell.bistratified( totalRGC, showPlots )
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) and cell
%   diameter (in mm) as a function of retinal eccentricity (in degrees).
%   Requires as input the output of the function cell.totalRGC.
%
% Examples:
%{
	totalRGC = cell.totalRGC();
	bistratified = cell.bistratified( totalRGC, true );
%}

% Handle plotting
if nargin==1
    showPlots = false;
end


%% Cell counts
% Data from Dacey 1993, Figure 13b:
%
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.
%
% Digitized data from Dacey's publications concerning morphology of
% bistratified ganglion cells allows us to calculate how much the volume of
% bistratified cells contributes to the overall RGC layer thickness from
% data concerning bistratified cells density varying with eccentricity
% across all meridians, and the cells' mean size.
%
% Dacey provided values for the temporal retina; we assume these values
% hold for all meridians

% Loop over the specified meridians
for mm = 1:length(totalRGC)

    % Set up this meridian model element
    bistratified(mm).label = totalRGC(mm).label;
    bistratified(mm).angle = totalRGC(mm).angle;

    % Data from Dacey 1993, used for all meridians
    supportMM = [0, .97, 1.96, 2.91, 3.92, 4.91, 5.92, 6.89, 7.84, 8.85, 9.86, 10.89, 11.9, 12.91, 13.84, 14.91];
    proportion = [0, .0135, .0168, .0202, .0241, .0284, .0324, .0364, .0403, .0447, .0485, .0538, .0573, .0603, .0641, .0662];

    % Convert retinal mm to visual degrees
    supportDeg = convert_mmRetina_to_degVisual(supportMM, bistratified(mm).angle);

    % Convert proportions to counts
    countsDegSq = proportion .* totalRGC(mm).countsDegSq(supportDeg);
    
    % Obtain a spline fit to the cell densities
    nonNan = ~isnan(countsDegSq);
    splineFit = fit(supportDeg(nonNan)', countsDegSq(nonNan)', 'cubicinterp');
    
    % Nan optic disc points and save the anonymous function
    bistratified(mm).countsDegSq = @(posDeg) ...        
        zeroOpticDiscPoints(splineFit(posDeg), posDeg, bistratified(mm).angle)';

    if showPlots
        if mm == 1
            figure
        end
        subplot(ceil(length(totalRGC)/2),ceil(length(totalRGC)/2),mm)
        plot(0:0.5:50,bistratified(mm).countsDegSq(0:0.5:50));
        hold on
        plot(supportDeg,countsDegSq,'*');
    end
    
end



end

