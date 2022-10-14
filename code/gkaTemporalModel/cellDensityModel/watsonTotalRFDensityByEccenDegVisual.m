function totalRFDensitySqDegVisual = watsonTotalRFDensityByEccenDegVisual(polarAngle)
% Implement Equation 8 of Watson 2014.
%
% Description:
%   This function returns the midget receptive field density (in units of
%   counts per square degree of the visual field), based upon the equation
%   and parameters provided for Equation 4 of Watson 2014.
%
%   NOTE: Watson's paper has the labels for the nasal and temporal retina
%   switched. In this routine we return the values that Watson labeled
%   "nasal" as the temporal retina, and vice-a-versa.
%
% Inputs:
%   supportPosDegVisual   - The positions (in degrees of visual angle) from
%                           the fovea at which to calculate the midget
%                           receptive field density
%   polarAngle            - The desired angle of the density function on
%                           the retinal field. (0=nasal; 90=superior;
%                           180=temporal; 270=inferior)
%
% Outputs:
%   totalRFDensitySqDegVisual - The density (receptive fields per visual
%                           square degree) of midget receptive fields at
%                           each of the positions
%

%% Check the input

if sum([0 90 180 270]==polarAngle) ~= 1
    error('The Watson equation for RF density is defined only for the cardinal meridia');
end

%% Obtain the parameters of the modeled receptive field density for this angle

% Taken from Table 1 of Watson
switch polarAngle
    case 0 % actual nasal (labeled temporal in Watson 2014)
        a = 0.9851;
        r2 = 1.058;
        re = 22.14;        
        
    case 90 % superior
        a = 0.9935;
        r2 = 1.035;
        re = 16.35;        
        
    case 180 % actual temporal (labeled nasal in Watson 2014)
        a = 0.9729;
        r2 = 1.084;
        re = 7.633;        
        
    case 270 % inferior
        a = 0.996;
        r2 = 0.9932;
        re = 12.13;        
        
end
        
% RGC fields (per deg visual ^2) at the foveal peak, same for every meridian
dgfZero = 33163.2;


% Note, the 2 * dcZero value is an expression for the number of midget
% retinal ganglion cells at the fovea, which is assumed to be 2 * dcZero.
% That is, it is assumed that the number of midget RGCs at the fovea is
% exactly equal to twice the number of cones.

watsonFunc = @(posDeg) dgfZero .* ...
    ( (a.*((1+posDeg./r2)).^-2) + (1-a).*exp(-1.*posDeg./re) );                                                    % The midget fraction at each eccentricity

% We do not zero out the optic disc points, as we will be integrating over
% the entire retina, for which the optic disc is just a tiny piece
%totalRFDensitySqDegVisual = @(posDeg) zeroOpticDiscPoints(watsonFunc(posDeg), posDeg, polarAngle);
totalRFDensitySqDegVisual = @(posDeg) watsonFunc(posDeg);

end

