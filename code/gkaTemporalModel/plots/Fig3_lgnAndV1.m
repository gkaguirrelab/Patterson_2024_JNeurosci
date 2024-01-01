% Region-level analysis of LGN and V1 responses. Bootstrap across
% acquisitions and use the Watson model fit to obtain the peak amplitude
% and frequency for each stimulus direction. Make plots of these.


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

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
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

% Define some ROI sets
roiSet = {'LGN','V1'};
nROIs = length(roiSet);

%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

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

% Load the LGN ROI
tmpPath = fullfile(localDataDir,'retinoFiles','LGN_bilateral.dtseries.nii');
LGNROI = cifti_read(tmpPath); LGNROI = LGNROI.cdata;


% Loop through subjects
for ss = 1:nSubs

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over bootstraps
    for bb = 1:nBoots

        % Get a sampling (with replacement) of the 12 acquisitions
        bootIdx = datasample(1:nAcqs,nAcqs);

        % Define some variables for parpool happiness
        nGood = []; peakFreq = []; peakAmp = []; goodIdx = []; Yset = [];

        % Loop over the ROIs
        for rr = 1:length(roiSet)

            switch roiSet{rr}
                case 'LGN'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (LGNROI == 1)));
                case 'V1'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 1)));
            end
            nGood(rr) = length(goodIdx);

            % Loop over the stimuli
            for whichStim = 1:nStims

                % Loop through the stimuli and frequencies and obtain the data
                rawVals = [];
                adjustedVals = [];

                for ff = 1:length(allFreqs)
                    subString = sprintf(['f%dHz_' stimulusDirections{whichStim}],allFreqs(ff));
                    idx = find(contains(stimLabels,subString));
                    rawVals(ff,:,:) = results.params(goodIdx,idx);
                end

                % Adjust the values for the zero frequency
                for ff = 2:length(allFreqs)
                    adjustedVals(ff-1,:,:) = rawVals(ff,:,:)-rawVals(1,:,:);
                end

                % Take the mean across voxels
                adjustedVals = squeeze(mean(adjustedVals,2));

                % Get the mean across (bootstrap resampled) acquisitions
                W = 1./std(adjustedVals(:,bootIdx),0,2)';
                Y = mean(adjustedVals(:,bootIdx),2)';

                % Perform the fit
                [p,~,~,yFitInterp] = fitWatsonModel(Y,W,studiedFreqs);

                % Determine the peak frequency in the log domain
                peakFreq(whichStim,rr) = mean(log10(interpFreqs(yFitInterp==max(yFitInterp))));

                % Save the peak amplitude, which is given by the first
                % param value
                peakAmp(whichStim,rr) = p(1);

                % Save the Y values
                Yset(whichStim,rr,:) = Y;

            end % stimuli

        end % ROIs

        % Store the bootstrap result in a par cell variable
        par_peakFreq{bb} = peakFreq;
        par_peakAmp{bb} = peakAmp;
        par_Yset{bb} = Yset;

    end % Bootstraps

    % Reshape the par cell data into matrices and get mean and SEM
    peakFreq =  reshape(cell2mat(par_peakFreq),nStims,nROIs,nBoots);
    peakFreqIQR(ss,:,:) = 10.^iqr(peakFreq,3);
    peakFreqMedian(ss,:,:) = 10.^median(peakFreq,3);

    peakAmp =  reshape(cell2mat(par_peakAmp),nStims,nROIs,nBoots);
    peakAmpIQR(ss,:,:) = iqr(peakAmp,3);
    peakAmpMedian(ss,:,:) = median(peakAmp,3);

    bootY = reshape(cell2mat(par_Yset),nStims,nROIs,nBoots,nFreqs);
    YMedian(ss,:,:,:) = squeeze(median(bootY,3));
    YIQR(ss,:,:,:) = squeeze(iqr(bootY,3));

end


%% Plot the TTFs

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [1 5.5 8]; % shifts each ttf down so they can be presented tightly on the same figure

% Prepare the TTF figure
figHandleA = figure('Renderer','painters');
figuresize(460,600,'pt');
tiledlayout(1,4,'TileSpacing','tight','Padding','tight')

% Loop over subjects
for whichSub = 1:length(subjects)

    % Loop over ROIS
    for rr = 1:nROIs

        % Select the plot of the correct stimulus direction
        nexttile((whichSub-1)*nROIs+rr);

        % Loop over stimuli and plot
        for whichStim = 1:nStims

            % Assemble the data
            thisY = squeeze(YMedian(whichSub,whichStim,rr,:))';
            thisIQR = squeeze(YIQR(whichSub,whichStim,rr,:))';
            lowY = thisY - thisIQR/2;
            highY = thisY + thisIQR/2;

            % Get the Watson fit
            [~,~,~,yFitInterp] = fitWatsonModel(thisY,1./thisIQR,studiedFreqs);

            % Add reference lines
            if whichStim == 2
                plot(log10([1 1]),[0 2]-shift_ttf(stimOrder(whichStim)),'-k');
            else
                plot(log10([1 1]),[0 4]-shift_ttf(stimOrder(whichStim)),'-k');
            end
            hold on
            plot(log10([1 100]),[0 0]-shift_ttf(stimOrder(whichStim)),'-k');

            % Add a patch for the error
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ lowY-shift_ttf(stimOrder(whichStim)), fliplr(highY)-shift_ttf(stimOrder(whichStim)) ],...
                plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha);

            % Add the model fit
            plot(log10(interpFreqs),yFitInterp-shift_ttf(stimOrder(whichStim)),...
                ['-' lineColor{stimOrder(whichStim)}],...
                'LineWidth',2);

            % Add the data symbols, using reversed markers for values below
            % zero
            idx = thisY > 0;
            plot(log10(studiedFreqs(idx)),thisY(idx)-shift_ttf(stimOrder(whichStim)),...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);
            idx = thisY < 0;
            plot(log10(studiedFreqs(idx)),thisY(idx)-shift_ttf(stimOrder(whichStim)),...
                'o','MarkerFaceColor','w',...
                'MarkerSize',6,'MarkerEdgeColor',lineColor{stimOrder(whichStim)},'LineWidth',1);

            % Clean up
            xlim(log10([0.5 150]))
            ylim([-8.5 5]);
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
            if whichSub~=1 || rr~=1
                a.XAxis.Visible = 'off';
            end
            title([subjects{whichSub} ' - ' roiSet{rr}]);

        end

        plot(log10([16 16]),[-8 5],'--k');

    end

end

% Save the plot
plotNamesPDF = 'Fig3a_lgnAndV1_withWatsonModel.pdf';
saveas(figHandleA,fullfile(savePath,plotNamesPDF));





%% Plot the peak Amplitude and Frequency


% Prepare the peak freq figure
figHandleB = figure('Renderer','painters');
figuresize(460,600,'pt');
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')
roiShift = 0.25;
% Loop over subjects
for whichSub = 1:length(subjects)

    % Select the plot of the correct stimulus direction
    nexttile();

    for whichStim = 1:nStims

        val = squeeze(peakAmpMedian(whichSub,whichStim,1));
        valIQR = squeeze(peakAmpIQR(whichSub,whichStim,1));

        % Loop over stimuli and plot
        plot([stimOrder(whichStim) stimOrder(whichStim)],[val-valIQR/2,val+valIQR/2],'-','Color',plotColor{stimOrder(whichStim)});
        hold on
        plot(stimOrder(whichStim),val,'o','Color',plotColor{stimOrder(whichStim)});

        val = squeeze(peakAmpMedian(whichSub,whichStim,2));
        valIQR = squeeze(peakAmpIQR(whichSub,whichStim,2));

        % Loop over stimuli and plot
        plot([stimOrder(whichStim)+roiShift stimOrder(whichStim)+roiShift],[val-valIQR/2,val+valIQR/2],'-','Color',plotColor{stimOrder(whichStim)});
        hold on
        plot(stimOrder(whichStim)+roiShift,val,'^','Color',plotColor{stimOrder(whichStim)});

    end

    ylim([0 6])
    ylabel('BOLD response [Pct change]');
    xlim([0.5 3.5]);
    a = gca;
    a.XTick = [1 1.25 2 2.25 3 3.25];
    a.XTickLabelRotation = 90;
    a.XTickLabel = {'lgn','v1','lgn','v1','lgn','v1'};

    title(subjects{whichSub});

end

% Loop over subjects
for whichSub = 1:length(subjects)

    % Select the plot of the correct stimulus direction
    nexttile();

    for whichStim = 1:nStims

        val = squeeze(peakFreqMedian(whichSub,whichStim,1));
        valIQR = squeeze(peakFreqIQR(whichSub,whichStim,1));

        % Loop over stimuli and plot
        plot([stimOrder(whichStim) stimOrder(whichStim)],[val-valIQR/2,val+valIQR/2],'-','Color',plotColor{stimOrder(whichStim)});
        hold on
        plot(stimOrder(whichStim),val,'o','Color',plotColor{stimOrder(whichStim)});

        val = squeeze(peakFreqMedian(whichSub,whichStim,2));
        valIQR = squeeze(peakFreqIQR(whichSub,whichStim,2));

        % Loop over stimuli and plot
        plot([stimOrder(whichStim)+roiShift stimOrder(whichStim)+roiShift],[val-valIQR/2,val+valIQR/2],'-','Color',plotColor{stimOrder(whichStim)});
        hold on
        plot(stimOrder(whichStim)+roiShift,val,'^','Color',plotColor{stimOrder(whichStim)});

    end

    ylim([0 30])
    ylabel('Peak frequency [Hz]');
    xlim([0.5 3.5]);

    a = gca;
    a.XTick = [1 1.25 2 2.25 3 3.25];
    a.XTickLabelRotation = 90;
    a.XTickLabel = {'lgn','v1','lgn','v1','lgn','v1'};

    title(subjects{whichSub});

end

% Save the plot
plotNamesPDF = 'Fig3b_lgnAndV1_withWatsonModel.pdf';
saveas(figHandleB,fullfile(savePath,plotNamesPDF));

% Report the amplification of the chromatic channels between LGN and V1,
% relative to the achromatic
for whichSub = 1:2
        denomer = squeeze(peakAmpMedian(whichSub,3,2))./squeeze(peakAmpMedian(whichSub,3,1));
    for whichStim = 1:2
        numer = squeeze(peakAmpMedian(whichSub,whichStim,2))./squeeze(peakAmpMedian(whichSub,whichStim,1));
        fprintf([subjects{whichSub} ' - ' stimulusDirections{whichStim} ': %2.2f \n'],numer/denomer);
    end
end
