% scriptCreatePlots

% Housekeeping
clear

% Properties of which model to plot
freqsForPlotting = logspace(0,2,50);

% Load the MRI data
mriData = loadMRIResponseData();

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% Load the model fitting results
prefs.verbose = false;
projectDir = tbLocateProject('mriSinaiAnalysis',prefs);
resultFilePath = fullfile(projectDir,'data','cstResultsBootstrap.mat');
load(resultFilePath,'results');

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the V1 cortical bins
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% The number of acquisitions obtained for each measurement. Might want this
% if we are going to do some boot-strapping
nAcqs = 12;

% Fixed features of the model
nCells = 3; nParams = 3;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [0 3 6 9 11 13]; % shifts each ttf down so they can be presented tightly on the same figure

% Loop over subjects
for whichSub = 1:length(subjects)


    % Get the median of the bootstrap parameters
    p = median(results.(subjects{whichSub}).p);

    % Get the y-values across bootstraps
    Y = [];
    for ii = 1:length(results.(subjects{whichSub}).Y)
        Y(:,:,:,ii) = results.(subjects{whichSub}).Y{ii};
    end
    yIQR = iqr(Y,4);
    Y = median(Y,4);
    yLow = Y - yIQR/2;
    yHi = Y + yIQR/2;
    yPlot = returnFitAcrossEccen(p,stimulusDirections,studiedEccentricites,freqsForPlotting);


    % Get the modeled response
    [response, rfsAtEcc] = returnFitAcrossEccen(p,stimulusDirections,studiedEccentricites,freqsForPlotting);

    % Prepare the figures
    figHandles = figure('Renderer','painters');
    figuresize(600,400,'pt');
    tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

    % Loop over stimuli and plot
    for whichStim = 1:nStims

        % Get the average response across eccentricity
        for ee=1:nEccs

            % Assemble the data
            v1YThisEcc = squeeze(Y(whichStim,ee,:))';
            v1YThisEcclow = squeeze(yLow(whichStim,ee,:))';
            v1YThisEcchigh = squeeze(yHi(whichStim,ee,:))';

            % Select the plot of the correct stimulus direction
            nexttile(stimOrder(whichStim));

            % Add the model fit
            % If this is LMS, show the separate midget and parasol components
            lineSpecs = {'-','--'};
            if whichStim == 3
                for cc=1:2
                    ttfComplex = double(subs(rfsAtEcc{ee}{3}(cc),freqsForPlotting));
                    vec = abs(ttfComplex);
                    plot(log10(freqsForPlotting),vec-shift_ttf(ee),...
                        lineSpecs{cc},'Color',[0.5 0.5 0.5],...
                        'LineWidth',2);
                    hold on
                end
            else
            plot(log10(freqsForPlotting),squeeze(response(whichStim,ee,:))-shift_ttf(ee),...
                ['-' lineColor{stimOrder(whichStim)}],...
                'LineWidth',2);
            hold on
            end


            % Add reference lines
            if ee==1 && whichStim == 3
                plot(log10([1 1]),[0 2],'-k');
            end
            plot(log10([1 2]),[0 0]-shift_ttf(ee),':k');
            plot(log10([50 100]),[0 0]-shift_ttf(ee),':k');

        end

    end

    % Clean up
    for ss=1:nStims
        for ee = 1:nEccs
            nexttile(ss);
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
            if ss>1
                a.XAxis.Visible = 'off';
            end
        end
    end

    % Save the plots
    plotNamesPDF = [subjects{whichSub} '_showMigetParasolModelFit.pdf' ];
    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end

