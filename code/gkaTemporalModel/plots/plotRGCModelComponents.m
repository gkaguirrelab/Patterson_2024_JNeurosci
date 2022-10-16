% plotRGCModelComponents

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract some info from the stored model structure
pRGC = rgcTemporalModel.p;
cfCone = rgcTemporalModel.cfCone;
coneDelay = rgcTemporalModel.coneDelay;
LMRatio = rgcTemporalModel.LMRatio;
pFitByEccen = rgcTemporalModel.pFitByEccen;
blockParamNames = rgcTemporalModel.meta.blockParamNames;
eccFields = rgcTemporalModel.meta.eccFields;
eccBins = rgcTemporalModel.meta.eccBins;
lbBlock = rgcTemporalModel.meta.lbBlock;
ubBlock = rgcTemporalModel.meta.ubBlock;

nEccBands = length(eccFields);
nBlockParams = size(pRGC,1);

% Plot parasol bipolar response across eccentricity
myFreqs = logspace(log10(0.5),log10(100),101);
myEccs = 1:5:86;
for ee = 1:length(myEccs)

    eccDeg = myEccs(ee);

    % Obtain the RGC model parameters for this eccentricity
    pBlockMidget = [];
    pBlockParasol = [];
    for ii = 1:7
        pBlockMidget(ii) = rgcTemporalModel.pFitByEccen{ii,1}(eccDeg);
        pBlockParasol(ii) = rgcTemporalModel.pFitByEccen{ii,2}(eccDeg);
    end

    % Derive the RGC models
    [rfRGCParasolLum, rfLMCone, rfParasolBipolar] = ...
        parseParamsParasol(pBlockParasol, cfCone, coneDelay);
    [rfRGCMidgetChrom, rfRGCMidgetLum, ~, rfMidgetBipolar] = ...
        parseParamsMidget(pBlockMidget, cfCone, coneDelay, LMRatio, eccDeg);

    % Construct the response grid
    rfLMConeGrid(ee,:) = abs(double(subs(rfLMCone,myFreqs)));
    rfParasolBipolarGrid(ee,:) = abs(double(subs(rfParasolBipolar,myFreqs)));
    rfMidgetBipolarGrid(ee,:) = abs(double(subs(rfMidgetBipolar,myFreqs)));
    rfRGCParasolLumGrid(ee,:) = abs(double(subs(rfRGCParasolLum,myFreqs)));
    rfRGCMidgetLumGrid(ee,:) = abs(double(subs(rfRGCMidgetLum,myFreqs)));
    rfRGCMidgetChromGrid(ee,:) = abs(double(subs(rfRGCMidgetChrom,myFreqs)));

end

rfLMConeGrid = (rfLMConeGrid./max(rfLMConeGrid,[],2)).*(max(rfLMConeGrid(:)));
figure
map = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';
colormap(map)
[X,Y] = meshgrid(log10(myFreqs),myEccs);
subplot(2,3,1);
s = mesh(X,Y,rfLMConeGrid);
s.FaceColor = 'interp';
s.FaceAlpha = 0.5;
subplot(2,3,2);
s = mesh(X,Y,rfParasolBipolarGrid);
s.FaceColor = 'interp';
s.FaceAlpha = 0.5;

subplot(2,3,3);
mesh(X,Y,rfMidgetBipolarGrid)
subplot(2,3,4);
mesh(X,Y,rfRGCParasolLumGrid)
subplot(2,3,5);
mesh(X,Y,rfRGCMidgetLumGrid)
subplot(2,3,6);
mesh(X,Y,rfRGCMidgetChromGrid)

figure
colormap(map)
contourf(X,Y,rfRGCMidgetChromGrid,15,'LineWidth',2)
