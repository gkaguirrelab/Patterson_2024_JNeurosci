% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','mriTemporalModel.mat');
load(loadPath,'mriTemporalModel');

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

% Plot each eccentricity band
eccDegs = [];
for ee = 1:nEccBands

    % Extract the midget parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,1));
    eccBin = eccBins{ee};
    eccField = eccFields{ee};

    % The eccProportion parameter is used to determine the precise
    % eccentricity location within the range provided by the eccBin
    eccProportion = pBlock(8);
    eccDegs(ee) = eccBin(1)+eccProportion*(range(eccBin));

    % Plot the midget temporal RFs
    figHandle = figure();

    % Get the midget temporal RFs for LminusM defined by these parameters
    [chromaticCenterWeight,chromaticSurroundWeight] = returnMidgetChromaticWeights(eccDegs(ee),LMRatio);
    rfRGC = returnRGCRF(pBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);
    plotRF(rfRGC,figHandle,'-r');

    rfRGC = returnRGCRF(pBlock,cfCone,coneDelay,1,1);
    plotRF(rfRGC,figHandle,'-k');

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
    saveas(gcf,fullfile('~/Desktop',plotName));

end

% Plot the parameters vs. eccentricity and obtain params x eccentricity
% Loop across cells
for cc=1:3
    figure

    % Loop across the 7 params that vary with eccentricty
    for ii=1:nBlockParams-1

        % The values for this param across eccentricity
        y = squeeze(pRGC(ii,:,cc));

        % Plot these values and the fit
        subplot(4,2,ii)
        plot(eccDegs,y,'ok')
        hold on
        eccDegsFit = 0:90;
        plot(eccDegsFit,pFitByEccen{ii,cc}(eccDegsFit),'-r')
        title(blockParamNames{ii})
        xlabel('Eccentricity [deg]'); ylabel('param value')
        ylim([lbBlock(ii) ubBlock(ii)]);

    end
    switch cc
        case 1
            sgtitle('Midget parameters')
        case 2
            sgtitle('Parasol parameters')
        case 3
            sgtitle('Bistratified parameters')
    end
end

% Figure of successive stages of the IRF



%% plot the model fits to the MRI data

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriTemporalModel.meta.studiedFreqs;
eccDegVals = mriTemporalModel.meta.eccDegVals;
subjects = mriTemporalModel.meta.subjects;
stimulusDirection = mriTemporalModel.meta.stimulusDirection;
plotColor = mriTemporalModel.meta.plotColor;
nFixed = mriTemporalModel.meta.nFixed;
nEcc = length(eccDegVals);
freqsForPlotting = logspace(0,2,50);

for whichStim = 1:length(stimulusDirection)
    for whichSub = 1:length(subjects)

        figure
        pMRI = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).p;
        v1Y = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Y;
        v1Y = reshape(v1Y,length(studiedFreqs),1,nEcc);
        v1Eccentricity = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Eccentricity;
        lgnFreqX = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnFreqX;
        lgnY = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnY;

        % Loop over eccentricities
        for ee=1:nEcc
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1Y(:,1,ee)),['o' plotColor{whichStim}]);
            hold on
            pMRIBlock = [pMRI(1:5) pMRI(nFixed+ee) pMRI(nFixed+nEcc+ee)];
            yFit = returnV1TTFForEcc(stimulusDirection{whichStim},rgcTemporalModel,eccDegVals(ee),pMRIBlock,freqsForPlotting);
            semilogx(freqsForPlotting,yFit,['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimulusDirection{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
            ylim([-1 7]);
        end

        % Add the LGN response
        lgnTTFFit = returnlgnTTF(stimulusDirection{whichStim},rgcTemporalModel,pMRI,freqsForPlotting,v1Eccentricity);
        subplot(2,4,8)
        semilogx(lgnFreqX,lgnY,['o' plotColor{whichStim}]);
        hold on
        semilogx(freqsForPlotting,lgnTTFFit,['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimulusDirection{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

        % Plot the surround suppression index vs. eccentricity
        subplot(2,4,4)
        plot(log10(eccDegVals),pMRI(:,nFixed+1:nFixed+nEcc),['*' plotColor{whichStim}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Save the plot
        plotName = [stimulusDirection{whichStim} '_' subjects{whichSub} '_ModelFit.pdf' ];
        saveas(gcf,fullfile('~/Desktop',plotName));
    end
end
