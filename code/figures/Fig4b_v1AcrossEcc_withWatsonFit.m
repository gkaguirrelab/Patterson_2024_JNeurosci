% scriptCreatePlots


%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;
nBoots = 250;

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(studiedFreqs);
interpFreqs = logspace(log10(1),log10(100),501);
nAcqs = 12;

% Define the eccentricity properties of the analysis
eccenDivs = [0 90./(2.^(5:-1:0))];
for ii=1:length(eccenDivs)-1
    eccenBins{ii}=[eccenDivs(ii),eccenDivs(ii+1)];
end
nEccs = length(eccenBins);

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [0 3 6 9 11 13]; % shifts each ttf down so they can be presented tightly on the same figure

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Loop over subjects
for ss = 1:length(subjects)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Load the retino maps for this subject
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    vArea = cifti_read(tmpPath); vArea = vArea.cdata;
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_eccen.dscalar.nii']);
    eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_angle.dscalar.nii']);
    polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;

    % Prepare the figure
    figHandles = figure('Renderer','painters');
    figuresize(600,400,'pt');
    tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

    % Loop over eccentricities
    for rr = 1:nEccs

        % Find the goodIdx
        eccenRange = eccenBins{rr};
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2))  ));

        % Loop through the stimuli
        for whichStim = 1:nStims

            % Loop through the frequencies and obtain the set of values
            rawVals = cell(1,nFreqs+1);
            for ff = 1:nFreqs+1
                subString = sprintf(['f%dHz_' stimulusDirections{whichStim}],allFreqs(ff));
                idx = find(contains(stimLabels,subString));
                rawVals{ff} = mean(results.params(goodIdx,idx),'omitnan');
            end

            % Adjust the values for the zero frequency
            adjustedVals = [];
            for ff = 2:nFreqs+1
                adjustedVals(:,ff-1) = rawVals{ff}-rawVals{1};
            end

            % Bootstrap over acquisitions to obtain the median data values
            % and the IQR
            bootY = zeros(nBoots,length(studiedFreqs));
            for bb = 1:nBoots
                bootY(bb,:) = mean(adjustedVals(datasample(1:nAcqs,nAcqs),:),1);
            end
            YMedian = median(bootY,1);
            YIQR = iqr(bootY,1);

            % A default p0 search point
            p0 = [1.5 5 1.1 1.5];

            % Special case a few p0 situations. The Watson model is quite
            % sensitive to the initial guess and prone to local minima in
            % fitting
            if (rr == 1) && (ss == 1) && (whichStim == 1)
                p0 = [3.5848    6.6386    0.8967    0.9982];
            end
            if (rr == 5) && (ss == 1) && (whichStim == 3)
                p0 = [2.5825    2.1325    2.6702    1.5639];
            end

            % Fit the model
            [pData(ss,whichStim,rr,:),fVal,~,yFitInterp] = fitWatsonModel(YMedian,1./YIQR,studiedFreqs,p0,interpFreqs);

            % Select the plot of the correct stimulus direction
            nexttile(stimOrder(whichStim));

            % Add a patch for the error
            Ylow = YMedian - YIQR/2;
            Yhi = YMedian + YIQR/2;
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ Ylow-shift_ttf(rr), fliplr(Yhi)-shift_ttf(rr) ],...
                plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha);
            hold on

            % Add the model fit
            plot(log10(interpFreqs),yFitInterp-shift_ttf(rr),...
                ['-' lineColor{stimOrder(whichStim)}],...
                'LineWidth',2);

            % Add the data symbols, using reversed markers for values below
            % zero
            idx = YMedian > 0;
            plot(log10(studiedFreqs(idx)),YMedian(idx)-shift_ttf(rr),...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);
            idx = YMedian < 0;
            plot(log10(studiedFreqs(idx)),YMedian(idx)-shift_ttf(rr),...
                'o','MarkerFaceColor','w',...
                'MarkerSize',6,'MarkerEdgeColor',lineColor{stimOrder(whichStim)},'LineWidth',1);

            % Add reference lines
            if rr==1 && whichStim == 3
                plot(log10([1 1]),[0 2],'-k');
            end
            plot(log10([1 2]),[0 0]-shift_ttf(rr),':k');
            plot(log10([50 100]),[0 0]-shift_ttf(rr),':k');

        end

    end

    % Clean up
    for whichStim=1:nStims
        nexttile(whichStim);
        xlim(log10([0.5 150]))
        ylim([-14 5])
        xlabel('Frequency [Hz]')
        a=gca;
        a.YTick = [0,2,4,6];
        a.YTickLabel = {'0','2','4','6'};
        a.XTick = log10([2,4,8,16,32,64]);
        a.XTickLabel = {'2','4','8','16','32','64'};
        a.XTickLabelRotation = 0;
        a.XMinorTick = 'off';
        a.YAxis.Visible = 'off';
        box off
        if whichStim>1
            a.XAxis.Visible = 'off';
        end
    end

    % Save the plots
    plotNamesPDF = [subjects{ss} '_v1Ecc_withWatsonModel.pdf' ];
    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end





