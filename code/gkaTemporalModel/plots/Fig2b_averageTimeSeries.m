% This script creates an illustrative figure of the average time-series
% data and model fit over all of area V1.

clear
close all

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% These variables define the subject names and stimulus directions
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
shortNames = {'gka','asb','cgp'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};

% Loop through the subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
    load(filePath,'results')

    % Properties of the stimuli
    stimLabels = results.model.opts{6};

    % Create a figure
    figHandle=portraitFigure();

    % Plot the average time series and model fit
    for ii=1:3
        subplot(4,1,ii);
        plot([1 672],[0 0],':k')
        hold on;
        plot(results.data.avgSignal(1+(ii-1)*672:672+(ii-1)*672),'.','Color',[0.5 0.5 0.5]);
        plot(results.data.avgModelFit(1+(ii-1)*672:672+(ii-1)*672),['-' plotColors{ii}],'LineWidth',2);
        ylim([-3 3]);
        title(directions{ii})
        ylabel('BOLD % change');
        xlim([1 336*2]);
        box off
        set(gca,'TickDir','out');
    end
    xlabel('time [seconds]')

    % Add the stimulus structure
    subplot(4,1,4)
    Xa = results.model.inputs{2}{1}(1:7,1:336).*(1:7)';
    Xb = results.model.inputs{2}{7}(49:55,1:336).*(1:7)';
    X = [nansum(Xa) nansum(Xb)];
    plot(X,'.k');
    yticks(1:8)
    yticklabels({'0 Hz','2 Hz','4 Hz','8 Hz','16 Hz','32 Hz','64 Hz'});

    set(gca,'XTick',[])
    set(gca,'TickDir','out');
    xlim([1 336*2]);
    box off
    set(gca,'XColor','none')

    % Save the figure
    set(figHandle,'color','none');
    fileName = fullfile(savePath,[subjectNames{ss} '_AvgV1_plots.pdf']);
    print(fileName,'-dpdf')

end



%% Local Functions
function f=portraitFigure()
f = figure();
set(f,...
    'PaperPosition',[0 0 11 8.5],...
    'PaperSize',[11 8.5000],...
    'PaperOrientation','portrait',...
    'Position',[543 336 791 611],...
    'OuterPosition',[543 336 791 690],...
    'InnerPosition',[543 336 791 611]);

end
