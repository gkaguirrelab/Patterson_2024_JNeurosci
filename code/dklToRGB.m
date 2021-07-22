function rgb = dklToRGB(dkl)
% Converts ±1 weights on the post-receptoral mechanism to a RGB display
%
% Syntax:
%  rgb = dklToRGB(dkl)
%
% Description:
%
%   This code is modified from the DKLDemo created by David Brainard and
%   distributed as part of the Psychtoolbox (and requires the Psychtoolbox
%   for execution).
%
% Inputs:
%   dkl                   - An nx3 matrix of weights, in the range of ±1,
%                           corresponding to the post-receptoral directions
%                           of red-green, blue-yellow, and luminance.
% Outputs:
%   rgb                   - An nx3 matrix of RGB color values, in the range
%                           of 0-1.
%
% Examples:
%{
    [X,Y,Z] = sphere(40);
    XVec = X(:); YVec = Y(:); ZVec = Z(:);
    dkl = [XVec, YVec, ZVec];
    rgb = dklToRGB(dkl);
    posQuadrant = logical((XVec<=1e-6) .* (YVec>=-1e-6) .* (ZVec>=0));
    rgb = reshape(rgb,41,41,3);
    figure
    surf(X,Y,Z,rgb,'EdgeColor','none','FaceColor','interp');
    axis equal
    view(-135,15);
%}

gamutHeadroom = 0.95;
lumBoost = 1;

% Load the Stockman-Sharpe cone fundamentals and a "standard" calibration
% file
load T_cones_ss2
load T_ss2000_Y2
S_cones = S_cones_ss2;
T_cones = T_cones_ss2;
T_Y = 683*T_ss2000_Y2;
S_Y = S_ss2000_Y2;
T_Y = SplineCmf(S_Y,T_Y,S_cones);
cal = LoadCalFile('PTB3TestCal');
calLMS = SetSensorColorSpace(cal,T_cones,S_cones);
calLMS = SetGammaMethod(calLMS,1);
calLum = SetSensorColorSpace(cal,T_Y,S_Y);

%% Define background.
% I go a bit below the monitor mid-point, allowing for more saturated
% display colors and a bit more gamut to work with.

bgLMS = PrimaryToSensor(calLMS,[0.4 0.4 0.4]');

%% Basic transformation matrices.  ComputeDKL_M() does the work.
% Get matrix that transforms between incremental
% cone coordinates and DKL coordinates 
% (Lum, RG, S).
M_ConeIncToDKL = ComputeDKL_M(bgLMS,T_cones,T_Y);
M_DKLToConeInc = inv(M_ConeIncToDKL);

%% Find incremental cone directions corresponding to DKL directions
lumConeInc = M_DKLToConeInc*[1 0 0]';
rgConeInc = M_DKLToConeInc*[0 1 0]';
sConeInc = M_DKLToConeInc*[0 0 1]';

% These directions are not scaled in an interesting way,
% need to scale them.  Here we'll find units so that 
% a unit excursion in the two directions brings us to
% the edge of the monitor gamut, with a little headroom.
bgPrimary = SensorToPrimary(calLMS,bgLMS);
lumPrimaryInc = SensorToPrimary(calLMS,lumConeInc+bgLMS)-bgPrimary;
rgPrimaryInc = SensorToPrimary(calLMS,rgConeInc+bgLMS)-bgPrimary;
sPrimaryInc = SensorToPrimary(calLMS,sConeInc+bgLMS)-bgPrimary;

lumScale = MaximizeGamutContrast(lumPrimaryInc,bgPrimary);
rgScale = MaximizeGamutContrast(rgPrimaryInc,bgPrimary);
sScale = MaximizeGamutContrast(sPrimaryInc,bgPrimary);

lumConeInc = gamutHeadroom*lumScale*lumConeInc*lumBoost;
rgConeInc = gamutHeadroom*rgScale*rgConeInc;
sConeInc = gamutHeadroom*sScale*sConeInc;


%% Perform the conversion

% Separate the three chanels from the dkl variable
XVec = dkl(:,1)'; YVec = dkl(:,2)'; ZVec = dkl(:,3)';
imageLMS = bgLMS*ones(size(XVec))+rgConeInc*XVec+sConeInc*YVec+lumConeInc*ZVec;
[rgb,badIndex] = SensorToSettings(calLMS,imageLMS);
bgRGB = SensorToSettings(calLMS,bgLMS);
rgb(:,find(badIndex == 1)) = bgRGB(:,ones(size(find(badIndex == 1))));
rgb = rgb';

end
