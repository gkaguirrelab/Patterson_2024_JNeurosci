function [fitParams, logisticFunc] = daceyMidgetFractionByEccenDegVisual(varargin)
% Calculate the midget fraction as a function of eccentricity
%
% Description:
%   This routine returns the parameters of a logistic fit to the
%   measurement of the fraction of retinal ganglion cells that are midget
%   RGCs. The calculation is based upon the values provided in Dacey (1993)
%   J Neurosci.
%
% Inputs:
%   none
%
% Optional key/value pairs:
%  'minRatio'             - The minimum value to use for the logistic fit
%                           to the Dacey plot.
%  'maxRatio'             - The maximum value to be used for the logistic 
%                           fit to the Dacey plot.
%  'logitFitStartPoint'   - Initial values to be used for the slope and 
%                           inflection point values of the logisitic fit.
%
% Outputs:
%   midgetFraction        - A 1 x p vector that contains the fraction of
%                           midget receptive fields at each eccentricity
%                           location.


%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('logitFitStartPoint',[7, 22, 0.45, 0.95],@isnumeric);
p.addParameter('showPlots',false,@islogical);

% parse
p.parse(varargin{:})


% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minRatio,maxRatio,x) minRatio+(maxRatio-minRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y');

% Data taken from Figure 19B of:
%	Dacey 1993 J Neurosci,
%   The Mosaic of Midget Ganglion Cells in the Human Retina
% digitzed using the web plot digitizer tool
daceyDataSupportPosMmRetina = [0.98174998	2.02854114	3.015966729	3.978015168	4.956959961	6.011644785	6.96599527	7.975927587	9.039745576	10.0028378	10.99287287	11.98068988	13.01939167	14.00851341	14.99998369];
daceyMidgetFraction = [91.54203702	97.96884938	95.9768409	94.73766615	82.97154041	68.49580038	58.38473457	52.33172959	48.381962	48.34575552	49.36116774	47.82027236	44.92408057	44.88689554	47.55638914];

% Convert the Dacey support vector from mm to degrees
daceyDataSupportPosDegVisual = convert_mmRetina_to_degVisual(daceyDataSupportPosMmRetina,0);

% Convert the Dacey percent values to proportion
daceyMidgetFraction = daceyMidgetFraction / 100;

% Fit a logistic function
logisticFit = fit(daceyDataSupportPosDegVisual',daceyMidgetFraction',logisticFunc, ...
    'StartPoint',p.Results.logitFitStartPoint );

fitParams = [logisticFit.slope, logisticFit.inflect, logisticFit.minRatio, logisticFit.maxRatio];

% Plot
if p.Results.showPlots
    figure
    plot(daceyDataSupportPosDegVisual,daceyMidgetFraction,'*r');
    hold on
    plot(0:90,logisticFit(0:90),'-k')
end

end