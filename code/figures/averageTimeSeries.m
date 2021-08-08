% Make some figures of the fit to the average V1 time series

%


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig X - time series fit');
mkdir(resultsSaveDir);

% Load the results file
fileName = fullfile(localSaveDir,'resultsFiles','HEROgka1_mtSinai_AvgV1_results.mat');
load(fileName,'results');


% Properties of the stimuli
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
stimLabels = results.model.opts{6};


figHandle=portraitFigure();

subplot(4,1,1)
plot(results.data.datats(1:336),'.','Color',[0.5 0.5 0.5]);
ylim([-3 3]);
xlim([1 336]);
ylabel('BOLD % change');
xlabel('time [seconds]')
box off
set(gca,'TickDir','out');
title('HEROgka1 LminusM Acq 01')

% Make a plot of the first LminusM acquisition temporal design
subplot(4,1,2)
imagesc(results.model.inputs{2}{1}(1:8,1:336).*(1:8)');
yticks(1:8)
yticklabels({'0 Hz','2 Hz','4 Hz','8 Hz','16 Hz','32 Hz','64 Hz','attention'});
cmap = zeros(9,3);
cmap(1,:) = [1 1 1];
cmap(2,:) = [0 0 0];
cmap(9,:) = [1 0 1];
for ii=1:6
    cmap(ii+2,:) = [0.75-(ii-1)*0.125, 0.75-(ii-1)*0.125, 1];
end

colormap(cmap);
set(gca,'XTick',[])
set(gca,'TickDir','out');
xlim([1 336]);
box off
set(gca,'XColor','none')

subplot(4,1,3)
hrf = results.data.hrf;
k=results.model.inputs{2}{1}(1:8,1:336).*results.params(56,1:8)';
for ii=1:8
    n=conv2( k(ii,:)',hrf);
    n=n(1:336);
    plot(n,'-','Color',cmap(ii+1,:));
    hold on
end
box off
set(gca,'XTick',[])
ylabel('BOLD % change');
ylim([0 6]);
xlim([1 336]);
set(gca,'TickDir','out');
set(gca,'XTick',[])
set(gca,'XColor','none')


subplot(4,1,4)
plot(results.data.datats(1:336),'.','Color',[0.5 0.5 0.5]);
hold on
plot(results.data.modelts(1:336),'-r','LineWidth',2);
ylim([-3 3]);
xlim([1 336]);
ylabel('BOLD % change');
xlabel('time [seconds]')
box off
set(gca,'TickDir','out');
set(gca,'XTick',[])
set(gca,'XColor','none')



% Save the figure
set(figHandle,'color','none');
fileName = fullfile(resultsSaveDir,'ModelFit_SingleAcquisition.pdf');
print(fileName,'-dpdf');






for ii=1:3
    subplot(3,1,ii);
    plot(results.data.datats(1+(ii-1)*4032:4032+(ii-1)*4032),'.','Color',[0.5 0.5 0.5]);
    hold on;
    plot(results.data.modelts(1+(ii-1)*4032:4032+(ii-1)*4032),'-r','LineWidth',2);
    ylim([-3 3]);
    xlim([0 4100]);
    title(directions{ii})
    ylabel('BOLD % change');
box off
set(gca,'TickDir','out');
end
xlabel('time [seconds]')

% Save the figure
set(figHandle,'color','none');
fileName = fullfile(resultsSaveDir,'ModelFit_AllAcquisitions.pdf');
print(fileName,'-dpdf')




figHandle=portraitFigure();

for ii=1:3
    subplot(3,1,ii);
    plot(results.data.avgSignal(1+(ii-1)*672:672+(ii-1)*672),'.','Color',[0.5 0.5 0.5]);
    hold on;
    plot(results.data.avgModelFit(1+(ii-1)*672:672+(ii-1)*672),'-r','LineWidth',2);
    ylim([-3 3]);
    title(directions{ii})
    ylabel('BOLD % change');
xlim([1 336*2]);
box off
set(gca,'TickDir','out');
end
xlabel('time [seconds]')

% Save the figure
set(figHandle,'color','none');
fileName = fullfile(resultsSaveDir,'ModelFit_AvgAcquisition.pdf');
print(fileName,'-dpdf')








figHandle=portraitFigure();

for dd = 1:3
    subplot(1,3,dd);
    for ff = 1:length(freqs)
        subString = sprintf(['f%dHz_' directions{dd}],freqs(ff));
        idx = find(contains(stimLabels,subString));
        vals{ff} = results.params(56,idx);
    end
    for ff = 2:length(freqs)
        semilogx(zeros(1,length(vals{ff}))+freqs(ff),vals{ff}-vals{1},'.','Color',[0.5 0.5 0.5]);
        hold on
    end
    meanVals = cellfun(@(x) mean(x-vals{1}),vals(2:end));
    semilogx(freqs(2:end),meanVals,'ob','MarkerSize',10,'MarkerFaceColor','b')
    title(directions{dd})
    ylabel('BOLD % change');
    xlabel('frequency [Hz]');
    semilogx([1 64],[0 0],':k','LineWidth',1)
    ylim([-1 4]);
    xlim([1 128])
    set(gca,'TickDir','out');
    box off
end


% Save the figure
set(figHandle,'color','none');
fileName = fullfile(resultsSaveDir,'ResponseByChannel.pdf');
print(fileName,'-dpdf')




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
