% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
[midgetData, parasolData] = loadRGCResponseData();

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

    % The eccProportion parameter is used to determine the precise eccentricity
    % location within the range provided by the eccBin
    eccProportion = pBlock(8);
    eccDegs(ee) = eccBin(1)+eccProportion*(range(eccBin));

    % Get the midget temporal RFs defined by these parameters
    [rfMidgetChrom, rfMidgetLum, rfLMCone] = ...
        parseParamsMidget(pBlock, cfCone, coneDelay, LMRatio, eccDegs(ee));

    % Plot the midget temporal RFs
    figHandle = figure();
    plotRF(rfMidgetLum,figHandle,'-k');
    plotRF(rfMidgetChrom,figHandle,'-r');
    plotRF(rfLMCone,figHandle,'-g',3);
    subplot(3,1,1);
    loglog(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.g,'*r');
    loglog(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.g,'*k');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));
    subplot(3,1,2);
    semilogx(midgetData.(eccField).chromatic.f,midgetData.(eccField).chromatic.p,'*r');
    semilogx(midgetData.(eccField).luminance.f,midgetData.(eccField).luminance.p,'*k');

    % Save the plot
    plotName = ['midgetTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile('~/Desktop',plotName));

    % Extract the parasol parameter values for this eccentricity band
    pBlock = squeeze(pRGC(:,ee,2));

    % Get the parasol temporal RF defined by these parameters
    [rfParasaolLum, rfLMCone] = ...
        parseParamsParasol(pBlock, cfCone, coneDelay);

    % Plot the parasol temporal RF
    figHandle = figure();
    plotRF(rfParasaolLum,figHandle,'-b');
    plotRF(rfLMCone,figHandle,'-g',3);
    subplot(3,1,1);
    hold on
    loglog(parasolData.(eccField).luminance.f,parasolData.(eccField).luminance.g,'*b');
    title(sprintf('Eccentricity = %2.1f',eccDegs(ee)));

    plotName = ['parasolTemporalRF_' num2str(eccDegs(ee),2) '_ModelFit.pdf' ];
    saveas(gcf,fullfile('~/Desktop',plotName));

end

% Plot the parameters vs. eccentricity and obtain params x eccentricity
% Loop across cells
for cc=1:2
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
    end
end

% Figure of successive stages of the IRF



%% plot the model fits to the MRI data

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriTemporalModel.meta.studiedFreqs;
eccDegVals = mriTemporalModel.meta.eccDegVals;
subjects = mriTemporalModel.meta.subjects;
stimuli = mriTemporalModel.meta.stimuli;
plotColor = mriTemporalModel.meta.plotColor;
nFixed = mriTemporalModel.meta.nFixed;
nEcc = length(eccDegVals);
freqsForPlotting = logspace(0,2,50);

for whichStim = 1:length(stimuli)
    for whichSub = 1:length(subjects)

        figure
        pMRI = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).p;
        v1Y = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1Y;
        v1Y = reshape(v1Y,length(studiedFreqs),1,nEcc);
        v1Eccentricity = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.v1Eccentricity;
        lgnFreqX = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.lgnFreqX;
        lgnY = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).data.lgnY;
        modelType = mriTemporalModel.(stimuli{whichStim}).(subjects{whichSub}).modelType;

        % Loop over eccentricities
        for ee=1:nEcc
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1Y(:,1,ee)),['o' plotColor{whichStim}]);
            hold on
            pBlock = [pMRI(1) pMRI(3:5) pMRI(nFixed+ee) pMRI(nFixed+nEcc+ee)];
            switch modelType
                case 'chromatic'
                    yFit = returnV1ChromEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
                case 'luminance'
                    yFit = returnV1LumEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
            end
            semilogx(freqsForPlotting,yFit,['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimuli{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
            ylim([-1 7]);
        end

        % Add the LGN response
        switch modelType
            case 'chromatic'
                lgnTTFFit = returnlgnLumTTFFit(pMRI,freqsForPlotting,v1Eccentricity);
            case 'luminance'
                lgnTTFFit = returnlgnLumTTFFit(pMRI,freqsForPlotting,v1Eccentricity);
        end
        subplot(2,4,8)
        semilogx(lgnFreqX,lgnY,['o' plotColor{whichStim}]);
        hold on
        semilogx(freqsForPlotting,lgnTTFFit,['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimuli{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

        % Plot the surround suppression index vs. eccentricity
        subplot(2,4,4)
        plot(log10(eccDegVals),pMRI(:,nFixed+1:nFixed+nEcc),['*' plotColor{whichStim}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Save the plot
        plotName = [stimuli{whichStim} '_' subjects{whichSub} '_ModelFit.pdf' ];
        saveas(gcf,fullfile('~/Desktop',plotName));
    end
end
