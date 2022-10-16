% plotRGCModel

% Housekeeping
clear
close all

% Load the empirical RGC data
[midgetData, parasolData] = loadRGCResponseData();

% Load the RGC model parameters
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','temporalModelResults','temporalModel.mat');
load(loadPath,'temporalModel');

% Extract some info from the stored model structure
p = temporalModel.p;
cfCone = temporalModel.cfCone;
coneDelay = temporalModel.coneDelay;
LMRatio = temporalModel.LMRatio;
pFitByEccen = temporalModel.pFitByEccen;
blockParamNames = temporalModel.meta.blockParamNames;
eccFields = temporalModel.meta.eccFields;
eccBins = temporalModel.meta.eccBins;
lbBlock = temporalModel.meta.lbBlock;
ubBlock = temporalModel.meta.ubBlock;

nEccBands = length(eccFields);
nBlockParams = size(p,1);

% Plot each eccentricity band
eccDegs = [];
for ee = 1:nEccBands

    % Extract the midget parameter values for this eccentricity band
    pBlock = squeeze(p(:,ee,1));
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

    % Extract the parasol parameter values for this eccentricity band
    pBlock = squeeze(p(:,ee,2));

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

end

% Plot the parameters vs. eccentricity and obtain params x eccentricity

% Loop across cells
for cc=1:2
    figure

    % Loop across the 7 params that vary with eccentricty
    for ii=1:nBlockParams-1

        % The values for this param across eccentricity
        y = squeeze(p(ii,:,cc));

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


