clear
close all

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
stimAlphas = [0.05 0.05 0.1];
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Visual area labels
roiLabels = {'V1','V2','V3','hV4','VO1','VO2','V3a','V3b','LO1','LO2','IPS0','IPS1'};

% Save a template map variable so we can create new maps below
templateImage = cifti_read(tmpPath);

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;


    % Prepare the figures
    figHandle = figure('Renderer','painters');
    figuresize(200,400,'pt');
    tiledlayout(2,1,'TileSpacing','tight','Padding','tight')

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
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
    load(filePath,'fitResults')


    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        posRespIdx = zeros(size(results.R2));
        posRespIdx(fitResults.eccDeg > 0) = cellfun(@(x) sum(x(whichStim,:))>1, fitResults.Y(fitResults.eccDeg > 0));

        % Get the median peak freq by area
        for vv = 1:12
            goodIdx = find(logical( (results.R2 > r2Thresh) .* posRespIdx .* (vArea == vv)  ));
            peakFreq = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));
            peakFreq(peakFreq == 1) = nan;
            peakFreq(peakFreq > 40) = nan;
            freqByROI(ss,whichStim,vv) = median(peakFreq,'omitmissing');           
        end

        % Identify those voxels with a positive response and an overall R2
        % of greater than the threshold
        goodIdx = find(logical( (results.R2 > r2Thresh) .* posRespIdx  ));
        nGood = length(goodIdx);

        % Extract the peak frequency for these vertices
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Filter the mis-fit vertices in which the max or minimum frequency
        % was identified
        peakFreq(peakFreq == 1) = nan;
        peakFreq(peakFreq > 40) = nan;

        % Update the goodIdx
        goodIdx = find(~isnan(peakFreq));

        % save a peakFreq map
        newMap = templateImage;
        newMap.cdata = single(zeros(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = single(peakFreq(goodIdx));
        newMap = ciftiMakePseudoHemi(newMap);
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
        cifti_write(newMap, fileOut);

        % Plot freq vs eccentricity for V1
        nexttile(ss);
        goodIdx = find(logical( (results.R2 > r2Thresh) .* posRespIdx .* ~isnan(peakFreq) .* (vArea == 1)  ));
        x = log10(fitResults.eccDeg(goodIdx));
        x(x<0)=x(x<0)/10;
        [x, sortedIdx] = sort(x);
        v = peakFreq(goodIdx);
        v = v(sortedIdx);
        pHandle=scatter(x,v,10,'w','o',...
            'MarkerEdgeColor','none','MarkerFaceColor',stimPlotColors{whichStim},...
            'MarkerFaceAlpha',stimAlphas(whichStim));

        hold on
        x(x<0)=0;
        % bin x and get median v in each bin
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
        a = gca();
        xTickVals = [1,2.5,5,10,20,40,80];
        xTickLabels = {'<1','2.5','5','10','20','40','80'};
        a.XTick = log10(xTickVals);
        a.XTickLabels = xTickLabels;
        xlim([-0.2 2.0]);
        xlabel('Eccentricity [deg]');
        ylabel('Peak frequency [Hz]');
        box off
    end

end

plotNamesPDF = 'freqByVertexEccen.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

