
clear
close all

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

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

    % Define model components
    datats = results.data.datats';
    for ii = 1:36
        data{ii}=datats((ii-1)*336+1:ii*336);
    end
    stimulus = results.model.inputs{2};
    modelOpts = results.model.opts;
    tr = results.meta.tr;
    modelClass = 'mtSinaiAdapt';

    % Special case of handling the stimtime variable for cgp
    if ss==3
        modelOpts = [modelOpts(1:9) {modelOpts(10:end)}];
    end

    % Add the stimClassSet we wish to model
    stimClassSet = {'f2Hz_LMS','f4Hz_LMS','f8Hz_LMS','f16Hz_LMS','f32Hz_LMS','f64Hz_LMS'};
    modelOpts{end+1}='stimClassSet';
    modelOpts{end+1}=stimClassSet;

    % Call the forwardModel
    adaptResults = forwardModel(data,stimulus,tr,...
        'averageVoxels',true,...
        'verbose',false,...
        'modelClass',modelClass,...
        'modelOpts',modelOpts);

    % Save the results
    fileName = fullfile(savePath,[subjectNames{ss} '_AvgV1_adaptMtSinai.mat']);
    save(fileName,'adaptResults');

    % Plot the data and timeseries fits
    figHandle=portraitFigure();
    for whichStim=1:3
        subplot(4,1,whichStim);
        plot([1 672],[0 0],':k')
        hold on;
        plot(adaptResults.data.avgSignal(1+(whichStim-1)*672:672+(whichStim-1)*672),'.','Color',[0.5 0.5 0.5]);
        plot(results.data.avgModelFit(1+(whichStim-1)*672:672+(whichStim-1)*672),['-' plotColors{whichStim}],'LineWidth',2);
        ylim([-3 3]);
        title(directions{whichStim})
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

    drawnow

    % Save the figure
    set(figHandle,'color','none');
    fileName = fullfile(savePath,[subjectNames{ss} '_AdaptModel_AvgV1_plots.pdf']);
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

