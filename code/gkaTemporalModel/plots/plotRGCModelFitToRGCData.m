% scriptCreatePlots

% Housekeeping
clear
close all

% Where to save figures
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Load the empirical RGC data
rcgData = loadRGCResponseData();

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
eccDegs = rgcTemporalModel.meta.eccDegs;

lbBlock = rgcTemporalModel.meta.lbBlock;
ubBlock = rgcTemporalModel.meta.ubBlock;

nEccBands = length(eccFields);
nBlockParams = size(pRGC,1);

% Plot each eccentricity band
for ee = 1:nEccBands

    % Extract the midget parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,1));
    eccField = eccFields{ee};

    % Plot the midget temporal RFs
    figHandle = figure();

    % Get the midget temporal RFs for LminusM defined by these parameters
    [chromaticCenterWeight,chromaticSurroundWeight] = returnRGCChromaticWeights('midget','LminusM',eccDegs(ee),LMRatio);
    rfRGC = returnRGCRF(pBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);
    plotRF(rfRGC,figHandle,'-r');

    [rfRGC, ~, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);
    plotRF(rfRGC,figHandle,'-k');
    plotRF(rfCone,figHandle,'-g',3);

    subplot(3,1,1);
    loglog(rcgData.midget.(eccField).LminusM.f,rcgData.midget.(eccField).LminusM.g,'*r');
    loglog(rcgData.midget.(eccField).LMS.f,rcgData.midget.(eccField).LMS.g,'*k');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    subplot(3,1,2);
    semilogx(rcgData.midget.(eccField).LminusM.f,rcgData.midget.(eccField).LminusM.p,'*r');
    semilogx(rcgData.midget.(eccField).LMS.f,rcgData.midget.(eccField).LMS.p,'*k');

    % Save the plot
    plotName = ['midgetTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile('~/Desktop',plotName));

    % Extract the parasol parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,2));

    % Get the parasol temporal RF defined by these parameters
    [rfRGC, ~, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);

    % Plot the parasol temporal RF
    figHandle = figure();
    plotRF(rfRGC,figHandle,'-k');
    plotRF(rfCone,figHandle,'-g',3);
    subplot(3,1,1);
    hold on
    loglog(rcgData.parasol.(eccField).LMS.f,rcgData.parasol.(eccField).LMS.g,'*k');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));

    plotName = ['parasolTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile('~/Desktop',plotName));

    % Extract the bistratified parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,3));

    % Get the parasol temporal RF defined by these parameters
    [rfRGC, ~, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);

    % Plot the bistratified temporal RF
    figHandle = figure();
    plotRF(rfRGC,figHandle,'-b');
    plotRF(rfCone,figHandle,'-g',3);
    subplot(3,1,1);
    hold on
    loglog(rcgData.bistratified.(eccField).S.f,rcgData.bistratified.(eccField).S.g,'*b');
    subplot(3,1,2);
    hold on
    semilogx(rcgData.bistratified.(eccField).S.f,rcgData.bistratified.(eccField).S.p,'*b');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));

    plotName = ['bistratifiedTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

end
