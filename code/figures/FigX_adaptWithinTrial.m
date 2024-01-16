clear
close all

% Place to save figures and load the results of the adapt analysis
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
stimClassSetLabels = {'chromatic','achromatic'};
nSubs = length(subjects);
plotColor={[0.85 0.55 0.85],[0.75 0.75 0.75]};
lineColor={'k','m'};
subMarkers = {'^','square','o'};
subLines = {'-','--',':'};
dirPlotShift = 0.1;
studiedFreqs = [2,4,8,16,32,64];

figure

% Loop through subjects
for whichSet = 1:length(stimClassSetLabels)
    for ss = 1:length(subjectNames)

        % Load the results file for this subject
        fileName = fullfile(savePath,[subjectNames{ss} '_' stimClassSetLabels{whichSet} '_AvgV1_adaptMtSinai.mat']);
        load(fileName,'adaptResults')

        % Get the adaptation params
        pp(ss,:) = adaptResults.params(end-8:end-3);

        plot(log10(studiedFreqs),pp(ss,:),subMarkers{ss},...
            'MarkerEdgeColor','none','MarkerFaceColor',plotColor{whichSet},...
            'MarkerSize',10)
        hold on

    end

    % Add median (across subject) values
    y = median(pp);
    pH(whichSet) = plot(log10(studiedFreqs),y,'-','Color',plotColor{whichSet},'LineWidth',2);
    plot(log10(studiedFreqs),y,'o',...
        'MarkerFaceColor','none','MarkerEdgeColor',plotColor{whichSet},...
        'MarkerSize',20,'LineWidth',2)
end

% Clean up
xlabel('Frequency [Hz]')
a=gca;
a.XTick = log10([2,4,8,16,32,64]);
a.XTickLabel = {'2','4','8','16','32','64'};
a.XTickLabelRotation = 0;
a.XMinorTick = 'off';
ylabel('Time constant of within-trial adaptation')
xlabel('Stimulus frequency [Hz]')



%
% % Plot the data and timeseries fits
% figHandle=portraitFigure();
% for whichStim=1:3
%     subplot(4,1,whichStim);
%     plot([1 672],[0 0],':k')
%     hold on;
%     plot(adaptResults.data.avgSignal(1+(whichStim-1)*672:672+(whichStim-1)*672),'.','Color',[0.5 0.5 0.5]);
%     plot(results.data.avgModelFit(1+(whichStim-1)*672:672+(whichStim-1)*672),['-' plotColors{whichStim}],'LineWidth',2);
%     ylim([-3 3]);
%     title(directions{whichStim})
%     ylabel('BOLD % change');
%     xlim([1 336*2]);
%     box off
%     set(gca,'TickDir','out');
% end
% xlabel('time [seconds]')
%
% % Add the stimulus structure
% subplot(4,1,4)
% Xa = results.model.inputs{2}{1}(1:7,1:336).*(1:7)';
% Xb = results.model.inputs{2}{7}(49:55,1:336).*(1:7)';
% X = [nansum(Xa) nansum(Xb)];
% plot(X,'.k');
% yticks(1:8)
% yticklabels({'0 Hz','2 Hz','4 Hz','8 Hz','16 Hz','32 Hz','64 Hz'});
%
% set(gca,'XTick',[])
% set(gca,'TickDir','out');
% xlim([1 336*2]);
% box off
% set(gca,'XColor','none')
%
% drawnow
%
% % Save the figure
% set(figHandle,'color','none');
% fileName = fullfile(savePath,[subjectNames{ss} '_LMSAdaptModel_AvgV1_plots.pdf']);
% print(fileName,'-dpdf')

