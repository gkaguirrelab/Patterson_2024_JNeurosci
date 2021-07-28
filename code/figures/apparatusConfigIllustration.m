%% apparatusConfigIllustration
%
% Make an illustration of the eye wearing a +36 contact lens and viewing a
% the stimulus field

% Stimulus distance and radius in mm
stimDistance = 6;
stimRadius = 6;

% Strength of the contact lens
lensDiopters = 36;

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig 1 - experimental design');
mkdir(resultsSaveDir);

% Create a figure
figHandle = figure();

% Render the eye
sceneGeometry=createSceneGeometry();
surfaceSet = sceneGeometry.refraction.retinaToCamera;
surfaceSet = addIris(surfaceSet, 2.1, 'green');
opticalSystem = surfaceSet.opticalSystem;
numRows = size(opticalSystem,1);
opticalSystem = opticalSystem(sum(isnan(opticalSystem),2)~=size(opticalSystem,2),:);
notLens = ~contains(surfaceSet.surfaceLabels,'lens');
opticalSystem = opticalSystem(notLens,:);
surfaceSet.opticalSystem = opticalSystem;
surfaceSet.surfaceColors = surfaceSet.surfaceColors(notLens);
surfaceSet.surfaceLabels = surfaceSet.surfaceLabels(notLens);
plotOpticalSystem(surfaceSet);

% Add a contact lens to the illustration
eye = sceneGeometry.eye;
opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'tears', eye.meta.spectralDomain ));
opticalSystemIn = [opticalSystem; eye.cornea.S(end,:) eye.cornea.side(end) eye.cornea.boundingBox(end,:) 1 1];
surfaceLabels = [eye.cornea.label(end); eye.cornea.label(end)];
surfaceColors = [eye.cornea.plot.color(end); eye.cornea.plot.color(end)];
surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
surfaceColors = [surfaceColors; {[0 0 1]}; {[0 0 1]}];
opticalSystemOut = addContactLens(opticalSystemIn,lensDiopters,'cornealRotation',eye.cornea.rotation,'contactLensViewAngle',35);
surfaceSet.opticalSystem = opticalSystemOut;
surfaceSet.surfaceLabels = surfaceLabels;
surfaceSet.surfaceColors = surfaceColors;
plotOpticalSystem(surfaceSet,'newFigure',false,'surfaceAlpha',0.5);

% Add the stimulus field
C = [stimDistance,0,0] ;   % center of circle 
R = stimRadius;    % Radius of circle 
theta=0:0.01:2*pi ;
x=C(1)+zeros(size(theta));
y=C(2)+R*sin(theta) ;
z = C(3)+R*cos(theta);
patch(x,y,z,[1 0 0],'FaceAlpha',0.5,'EdgeColor','k','LineWidth',2);
hold on

% Clean up the axes
view(30,8)
box off
axis off

% Calculate the magnification produced by the contact lens
magEffect = calcAngularMagnification(eye,'contactLens',lensDiopters,'contactLensViewAngle',35);

% Report the radial width (in degrees) of the stimulus field
width = atand(stimRadius/stimDistance);
title(sprintf('radial stimulus width = %2.1f°, x%2.2f',width,magEffect));

% Save the figure
set(figHandle,'color','none');
fileName = fullfile(resultsSaveDir,'EyeAndStimulusField.png');
print(fileName,'-dpng','-r1200','-painters');

