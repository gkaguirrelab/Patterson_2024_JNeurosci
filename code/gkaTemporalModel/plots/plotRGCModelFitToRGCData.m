% scriptCreatePlots

% Housekeeping
clear
%close all

% Where to save figures
savePath = fullfile('~','Desktop','VSS 2023');

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
rgcTemporalModel = fitRGCFResponse();

% Extract some info from the stored model structure
pRGC = rgcTemporalModel.p;
cfCone = rgcTemporalModel.cfCone;
coneDelay = rgcTemporalModel.coneDelay;
LMRatio = rgcTemporalModel.LMRatio;
pFitByEccen = rgcTemporalModel.pFitByEccen;
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
    figuresize(200, 400,'pt');

    % Get the midget RFs
    [chromaticCenterWeight,chromaticSurroundWeight] = returnRGCChromaticWeights('midget','LminusM',eccDegs(ee),LMRatio);
    rfRGC = returnRGCRF(pBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);
    plotRF(rfRGC,figHandle,'-r');
    subplot(3,1,1);
    ylim([0 10]);

    [rfRGC, rfBipolar, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);
    plotRF(rfRGC,figHandle,'-k');
    subplot(3,1,1);
    ylim([0 10]);
    box off

    subplot(3,1,1);
    loglog(rcgData.midget.(eccField).LminusM.f,rcgData.midget.(eccField).LminusM.g,'or');
    loglog(rcgData.midget.(eccField).LMS.f,rcgData.midget.(eccField).LMS.g,'ok');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    ylim([0 10]);
    box off

    subplot(3,1,2);
    semilogx(rcgData.midget.(eccField).LminusM.f,rcgData.midget.(eccField).LminusM.p,'or');
    semilogx(rcgData.midget.(eccField).LMS.f,rcgData.midget.(eccField).LMS.p,'ok');
    box off

    % Save the plot
    plotName = ['midgetTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

    % Extract the parasol parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,2));

    % Get the parasol temporal RF defined by these parameters
    [rfRGC, ~, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);

    % Plot the parasol temporal RF
    figHandle = figure();
    figuresize(200, 400,'pt');

    plotRF(rfRGC,figHandle,'-k');
    box off

    subplot(3,1,1);
    hold on
    loglog(rcgData.parasol.(eccField).LMS.f,rcgData.parasol.(eccField).LMS.g,'ok');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    ylim([0 10]);
    box off
    subplot(3,1,2);
    box off

    plotName = ['parasolTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

    % Extract the bistratified parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,3));

    % Get the parasol temporal RF defined by these parameters
    [rfRGC, ~, rfCone] = returnRGCRF(pBlock,cfCone,coneDelay,1,1);

    % Plot the bistratified temporal RF
    figHandle = figure();
    figuresize(200, 400,'pt');

    plotRF(rfRGC,figHandle,'-b');
    box off

    subplot(3,1,1);
    hold on
    loglog(rcgData.bistratified.(eccField).S.f,rcgData.bistratified.(eccField).S.g,'ob');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    ylim([0 10]);
    box off

    plotName = ['bistratifiedTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

end
