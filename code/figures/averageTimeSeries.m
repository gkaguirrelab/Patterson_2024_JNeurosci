% This script creates an illustrative figure of the average time-series
% data and model fit over all of area V1.


% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};
analysisIDs = { '61d32407d1304b39ec5427dc','61d3242146c9dab1751bd55f' };


% Create a flywheel object. You need to set your flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig 2 - time series analysis');
mkdir(resultsSaveDir);

% Loop through the subjects
for ss = 1:length(subjectNames)

    % Download and then load the results file
    fileName = [subjectNames{ss} '_mtSinai_results.mat'];
    tmpPath = fullfile(resultsSaveDir,[subjectNames{ss} '_mtSinai_AvgV1_results.mat']);
    fw.downloadOutputFromAnalysis(analysisIDs{ss},fileName,tmpPath);
    load(tmpPath,'results');

    % Properties of the stimuli
    stimLabels = results.model.opts{6};

    % Create a figure
    figHandle=portraitFigure();

    % Plot the average time series and model fit
    for ii=1:3
        subplot(4,1,ii);
        plot(results.data.avgSignal(1+(ii-1)*672:672+(ii-1)*672),'.','Color',[0.5 0.5 0.5]);
        hold on;
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
    Xa = results.model.inputs{2}{1}(1:7,1:336).*(7:-1:1)';
    Xb = results.model.inputs{2}{7}(49:55,1:336).*(7:-1:1)';
    X = flipud([Xa Xb]);
  imagesc(X);
    yticks(1:8)
    yticklabels({'64 Hz','32 Hz','16 Hz','8 Hz','4 Hz','2 Hz','0 Hz'});
    cmap = zeros(8,3);
    cmap(1,:) = [1 1 1];
    for ii=2:7
        cmap(9-ii,:) = [0.75-(ii-1)*0.125, 0.75-(ii-1)*0.125, 1];
    end

    colormap(cmap);
    set(gca,'XTick',[])
    set(gca,'TickDir','out');
    xlim([1 336*2]);
    box off
    set(gca,'XColor','none')

    % Save the figure
    set(figHandle,'color','none');
    fileName = fullfile(resultsSaveDir,[subjectNames{ss} '_AvgV1_plots.pdf']);
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
