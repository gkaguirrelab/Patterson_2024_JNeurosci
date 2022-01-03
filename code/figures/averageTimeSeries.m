% Make some figures of the fit to the average V1 time series
%


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.25;

% This is the visual area and eccentricity range to grab. The visual areas
% are: V1 = 1, V2 = 2, V3 = 3, hV4/LO = [4 5], MT/MST = [8 9]
area = 1;

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig 2 - time series analysis');
mkdir(resultsSaveDir);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};

% Loop through the subjects
for ss = 1:length(subjectNames)

% Load the results file
fileName = fullfile(localSaveDir,'resultsFiles',[subjectNames{ss} '_mtSinai_AvgV1_results.mat']);
load(fileName,'results');

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

% Make a plot of the first LminusM acquisition temporal design
subplot(4,1,2)
imagesc(results.model.inputs{2}{1}(1:7,1:336).*(1:7)');
yticks(1:8)
yticklabels({'0 Hz','2 Hz','4 Hz','8 Hz','16 Hz','32 Hz','64 Hz'});
cmap = zeros(9,3);
cmap(1,:) = [1 1 1];
cmap(2,:) = [0 0 0];
%cmap(9,:) = [1 0 1];
for ii=1:6
    cmap(ii+2,:) = [0.75-(ii-1)*0.125, 0.75-(ii-1)*0.125, 1];
end

colormap(cmap);
set(gca,'XTick',[])
set(gca,'TickDir','out');
xlim([1 336]);
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
