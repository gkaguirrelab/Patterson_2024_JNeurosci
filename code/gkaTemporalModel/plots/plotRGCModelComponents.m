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
    plotData.rfLMConeGrid(ee,:) = abs(double(subs(rfLMCone,myFreqs)));
    plotData.rfParasolBipolarGrid(ee,:) = abs(double(subs(rfParasolBipolar,myFreqs)));
    plotData.rfMidgetBipolarGrid(ee,:) = abs(double(subs(rfMidgetBipolar,myFreqs)));
    plotData.rfRGCParasolLumGrid(ee,:) = abs(double(subs(rfRGCParasolLum,myFreqs)));
    plotData.rfRGCMidgetLumGrid(ee,:) = abs(double(subs(rfRGCMidgetLum,myFreqs)));
    plotData.rfRGCMidgetChromGrid(ee,:) = abs(double(subs(rfRGCMidgetChrom,myFreqs)));

end

plotItems = {'rfLMConeGrid','rfParasolBipolarGrid','rfMidgetBipolarGrid',...
    'rfRGCParasolLumGrid','rfRGCMidgetLumGrid','rfRGCMidgetChromGrid'};
[X,Y] = meshgrid(log10(myFreqs),myEccs);
zLimMax = [60,40,40,10,10,10];
map = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';

f = figure('Renderer','painters');
colormap(map)
for ii=1:length(plotItems)
    subplot(2,3,ii);
    plot3([log10(myFreqs(1)) log10(myFreqs(1))],[0 90],[0 0],'-b')
    hold on
    plot3([log10(myFreqs(1)) log10(myFreqs(1))],[0 90],[zLimMax(ii) zLimMax(ii)],'-b')
    plot3([log10(myFreqs(end)) log10(myFreqs(end))],[0 90],[zLimMax(ii) zLimMax(ii)],'-b')
    plot3([log10(myFreqs(1)) log10(myFreqs(1))],[0 0],[0 zLimMax(ii)],'-b')
    plot3([log10(myFreqs(1)) log10(myFreqs(1))],[90 90],[0 zLimMax(ii)],'-b')
    plot3([log10(myFreqs(1)) log10(myFreqs(end))],[0 0],[0 0],'-b')
    plot3([log10(myFreqs(1)) log10(myFreqs(end))],[0 0],[zLimMax(ii) zLimMax(ii)],'-b')
    plot3([log10(myFreqs(1)) log10(myFreqs(end))],[90 90],[zLimMax(ii) zLimMax(ii)],'-b')
    s = mesh(X,Y,plotData.(plotItems{ii}));
    plot3([log10(myFreqs(end)) log10(myFreqs(end))],[90 90],[0 zLimMax(ii)],'-b')
    s.FaceColor = 'interp';
    s.FaceAlpha = 0.5;
    s.EdgeAlpha = 0.75;
    view([145 20]);
    ylim([0 90]);
    zlim([0 zLimMax(ii)]);
    a = gca;
    a.Color = 'none';
    a.XTickLabel={'1','10','100'};
    a.YTick=[0 30 60 90];
    a.YTickLabelRotation = 0;
    axis square
    grid off
    box off
    title({plotItems{ii}})
end

saveas(f,'~/Desktop/untitled.pdf');

figure
colormap(map)
for ii=1:length(plotItems)
    subplot(2,3,ii);
contourf(X,Y,plotData.(plotItems{ii}),15,'LineWidth',0.5)
    axis square
end
