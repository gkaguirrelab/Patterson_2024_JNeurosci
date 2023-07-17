clear
close all

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
stimAlphas = [0.05 0.05 0.1];
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


% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(400,400,'pt');
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(savePath,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFit_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValSet(fitResults.eccDeg > 0) = cellfun(@(x) x(whichStim), fitResults.fVal(fitResults.eccDeg > 0));
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));

        % Extract the peak amplitude and frequency for these vertices
        peakAmp = nan(nVert,1);
        peakAmp(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakAmp(goodIdx));
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Plot Amplitude vs eccentricity for V1
        nexttile((ss-1)*2+1);
        x = log10(fitResults.eccDeg(goodIdx));
        x(x<0)=x(x<0)/10;
        [x, sortedIdx] = sort(x);
        v = peakAmp(goodIdx);
        v = v(sortedIdx);
        v(v>6) = 6 + (v(v>6)-6)/5;
        pHandle=scatter(x,v,10,'w','o',...
            'MarkerEdgeColor','none','MarkerFaceColor',stimPlotColors{whichStim},...
            'MarkerFaceAlpha',stimAlphas(whichStim));

        % Add a fit line through the median value in each of 30
        % eccentricity bins
        hold on
        x(x<0)=0;
        nBins = 30;
        [binIdx,edges] = discretize(x,nBins);
        binCenters = edges(1:end-1)+diff(edges)/2;
        for ii = 1:nBins
            binV(ii) = median(v(binIdx==ii));
        end
        % plot(binCenters,binV,'*','Color',stimPlotColors{whichStim});
        sp = spaps(binCenters,binV,-750);
        xq = binCenters(1):0.01:1.8;
        vq = fnval(sp,xq);
        plot(xq,vq,['-' stimPlotColors{whichStim}],'LineWidth',1.5)

        % Clean up
        ylim([0 7]);
        a = gca();
        a.YTick = [0 1 2 3 4 5 6];
        a.YTickLabels = {'0','1','2','3','4','5','>6'};
        xTickVals = [1,2,5,10,20,40,80];
        xTickLabels = {'<1','2','5','10','20','40','80'};
        a.XTick = log10(xTickVals);
        a.XTickLabels = xTickLabels;
        xlim([-0.2 2.0]);
        xlabel('Eccentricity [deg]');
        ylabel('Amplitude [BOLD Pct change]');
        box off

        title(subjects{ss});


        % Plot freq vs eccentricity for V1
        nexttile((ss-1)*2+2);
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));
        x = log10(fitResults.eccDeg(goodIdx));
        x(x<0)=x(x<0)/10;
        [x, sortedIdx] = sort(x);
        v = peakFreq(goodIdx);
        v = v(sortedIdx);
        v(v>40) = 40 + (v(v>40)-40)/10;
        pHandle=scatter(x,v,10,'w','o',...
            'MarkerEdgeColor','none','MarkerFaceColor',stimPlotColors{whichStim},...
            'MarkerFaceAlpha',stimAlphas(whichStim));

        % Add a fit line through the median value in each of 30
        % eccentricity bins
        hold on
        x(x<0)=0;
        nBins = 30;
        [binIdx,edges] = discretize(x,nBins);
        binCenters = edges(1:end-1)+diff(edges)/2;
        for ii = 1:nBins
            binV(ii) = median(v(binIdx==ii));
        end
        sp = spaps(binCenters,binV,-750);
        xq = binCenters(1):0.01:1.8;
        vq = fnval(sp,xq);
        plot(xq,vq,['-' stimPlotColors{whichStim}],'LineWidth',1.5)

        % Clean up
        ylim([0 45]);
        a = gca();
        a.YTick = [0 10 20 30 40];
        a.YTickLabels = {'0','10','20','30','>40'};
        xTickVals = [1,2,5,10,20,40,80];
        xTickLabels = {'<1','2','5','10','20','40','80'};
        a.XTick = log10(xTickVals);
        a.XTickLabels = xTickLabels;
        xlim([-0.2 2.0]);
        xlabel('Eccentricity [deg]');
        ylabel('Peak frequency [Hz]');
        box off

        title(subjects{ss});


    end

end

plotNamesPDF = 'ampAndFreqByVertexEccen.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

