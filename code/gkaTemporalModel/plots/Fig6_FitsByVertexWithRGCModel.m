clear
close all

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
subMarkers = {'^','square'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Save a template map variable so we can create new maps below
templateImage = cifti_read(tmpPath);

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.1;

% This is the threshold for the goodness of fit of the Watson model to the
% TTF in each vertex. We only display those voxels with this fVal or lower
fValThresh = 2;

% Define the eccentricity locations for the rgc model
modelEccentricities = [5 21 43];
nEccs = length(modelEccentricities);

% The number of bins for the fMRI data across eccentricity
nEccDataBins = 30;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];
interpFreqs = logspace(log10(1),log10(100),501);

% Create a p vector that produces no post-retinal modification
nCells = 3; nParams = 3;
pMRI = [1 repmat([200 1 0.15],1,nCells*nEccs)];

% Get the rgc modeled response
rgcResponse = returnCSTFitAcrossEccen(pMRI,stimulusDirections,modelEccentricities,interpFreqs);

% Extract the maximum amplitude and peak frequency across stimuli and
% eccentricites
for whichStim = 1:nStims
    for ee=1:nEccs
        Y = squeeze(rgcResponse(whichStim,ee,:));
        peakRGCAmp(whichStim,ee) = max(Y);
        idx = find(Y == max(Y));
        peakRGCFreq(whichStim,ee) = interpFreqs(idx(1));
    end
end

% Scale the amplitude of RGC response to the mean LMS response
peakRGCAmp = peakRGCAmp ./ mean(peakRGCAmp(3,:));

% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(400,200,'pt');
tiledlayout(1,2,'Padding','tight')

% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    ampNorm = [];

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(savePath,[subjectNames{ss} '_WatsonFit_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValIdx = find(cellfun(@(x) ~isempty(x), fitResults.fVal));
        fValSet(fValIdx) = cellfun(@(x) x(whichStim), fitResults.fVal(fValIdx));
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));

        % Extract the peak amplitude and frequency for these vertices
        peakAmp = nan(nVert,1);
        peakAmp(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakAmp(goodIdx));
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Plot Amplitude vs eccentricity for V1
        nexttile(1);
        x = log10(eccenMap(goodIdx));
        x(x<0)=x(x<0)/10;
        [x, sortedIdx] = sort(x);
        v = peakAmp(goodIdx);
        v = v(sortedIdx);
        v(v>6) = 6 + (v(v>6)-6)/5;

        % Add a fit line through the median value in each of 30
        % eccentricity bins
        hold on
        x(x<0)=0;
        [binIdx,edges] = discretize(x,nEccDataBins);
        binCenters = edges(1:end-1)+diff(edges)/2;
        for ii = 1:nEccDataBins
            binV(ii) = median(v(binIdx==ii));
        end
        % plot(binCenters,binV,'*','Color',stimPlotColors{whichStim});
        sp = spaps(binCenters,binV,-750);
        xq = binCenters(1):0.01:1.8;
        vq = fnval(sp,xq);

        xq = binCenters;
        vq = binV;

        % Normalize by RGC response to LMS
        if whichStim == 3
            ampNorm = vq;
        end
        vq = vq ./ ampNorm;

        scatter(xq,log10(vq),subMarkers{ss},stimPlotColors{whichStim},...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor',stimPlotColors{whichStim},...
            'MarkerFaceAlpha',0.2)
        hold on

        % Add the rgc model line. Only extend to 18 degrees for the
        % bistratified model
        if whichStim == 2
            plot(log10(modelEccentricities(1:2)),log10(peakRGCAmp(whichStim,1:2)),...
                '*-','Color',stimPlotColors{whichStim},'LineWidth',2);
        else
            plot(log10(modelEccentricities),log10(peakRGCAmp(whichStim,:)),...
                '*-','Color',stimPlotColors{whichStim},'LineWidth',2);
        end

        % Clean up
        a = gca();
        a.YTick = [-2 -1 0 1];
        a.YTickLabels = {'0.01','0.1','1','10'};
        xTickVals = [1,2,5,10,20,40,80];
        xTickLabels = {'1','2','5','10','20','40','80'};
        a.XTick = log10(xTickVals);
        a.XTickLabels = xTickLabels;
        xlim([-0.2 2.0]);
        xlabel('Eccentricity [deg]');
        ylabel('Relative amplitude [a.u.]');
        box off

        title('peak amplitude');


        % Plot freq vs eccentricity for V1
        nexttile(2);
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));
        x = log10(eccenMap(goodIdx));
        x(x<0)=x(x<0)/10;
        [x, sortedIdx] = sort(x);
        v = peakFreq(goodIdx);
        v = v(sortedIdx);
        v(v>40) = 40 + (v(v>40)-40)/10;

        % Add a fit line through the median value in each of 30
        % eccentricity bins
        hold on
        x(x<0)=0;
        nEccDataBins = 30;
        [binIdx,edges] = discretize(x,nEccDataBins);
        binCenters = edges(1:end-1)+diff(edges)/2;
        for ii = 1:nEccDataBins
            binV(ii) = median(v(binIdx==ii));
        end
        sp = spaps(binCenters,binV,-750);
        xq = binCenters(1):0.01:1.8;
        vq = fnval(sp,xq);

        xq = binCenters;
        vq = binV;

        scatter(xq,vq,subMarkers{ss},stimPlotColors{whichStim},...
            'MarkerEdgeColor','none',...
            'MarkerFaceColor',stimPlotColors{whichStim},...
            'MarkerFaceAlpha',0.2)

        % Add the rgc model line. Only extend to 18 degrees for the
        % bistratified model
        if whichStim == 2
            plot(log10(modelEccentricities(1:2)),peakRGCFreq(whichStim,1:2),...
                '*-','Color',stimPlotColors{whichStim},'LineWidth',2);
        else
            plot(log10(modelEccentricities),peakRGCFreq(whichStim,:),...
                '*-','Color',stimPlotColors{whichStim},'LineWidth',2);
        end

        % Clean up
        a = gca();
        xTickVals = [1,2,5,10,20,40,80];
        xTickLabels = {'1','2','5','10','20','40','80'};
        a.XTick = log10(xTickVals);
        a.XTickLabels = xTickLabels;
        xlim([-0.2 2.0]);
        xlabel('Eccentricity [deg]');
        ylabel('Peak frequency [Hz]');
        box off

        title('peak frequency');


    end

end

plotNamesPDF = 'Fig6_ampAndFreqVsRGCModel.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

