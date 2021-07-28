eye = modelEyeParameters('sphericalAmetropia',0);

plotOpticalSystem(eye);

opticalSystem = initializeOpticalSystem(returnRefractiveIndex( 'tears', eye.meta.spectralDomain ));

% Add the tear surface
opticalSystemIn = [opticalSystem; eye.cornea.S(end,:) eye.cornea.side(end) eye.cornea.boundingBox(end,:) 1 1];

% Assemble the labels
surfaceLabels = [eye.cornea.label(end); eye.cornea.label(end)];

% Assemble the surface plot colors
surfaceColors = [eye.cornea.plot.color(end); eye.cornea.plot.color(end)];

surfaceLabels = [surfaceLabels; {'contactLens'}; {'tearfilm'}];
surfaceColors = [surfaceColors; {[0 0 1]}; {[0 0 1]}];


lensDiopters = 36;
eyePower = calcOpticalPower(opticalSystemIn);
opticalSystemOut = addContactLens(opticalSystemIn,lensDiopters,'cornealRotation',eye.cornea.rotation,'contactLensViewAngle',35);

surfaceSet.opticalSystem = opticalSystemOut;
surfaceSet.surfaceLabels = surfaceLabels;
surfaceSet.surfaceColors = surfaceColors;

plotOpticalSystem(surfaceSet,'newFigure',false,'surfaceAlpha',0.5);

C = [12,0,0] ;   % center of circle 
R = 12. ;    % Radius of circle 
theta=0:0.01:2*pi ;
x=C(1)+zeros(size(theta));
y=C(2)+R*sin(theta) ;
z = C(3)+R*cos(theta);
patch(x,y,z,[0.5 1 0.5],'FaceAlpha',0.5)
hold on


calcAngularMagnification(eye,'contactLens',lensDiopters,'contactLensViewAngle',35)

