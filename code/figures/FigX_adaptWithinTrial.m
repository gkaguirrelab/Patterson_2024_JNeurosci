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
stimClassParamIdx = {[289:294],[295:300]};
nSubs = length(subjects);
plotColor={[0.85 0.55 0.85],[0.75 0.75 0.75]};
lineColor={'m','k'};
subMarkers = {'^','square','o'};
dirPlotShift = 0.03;
studiedFreqs = [2,4,8,16,32,64];

figHandle = figure();
figuresize(400,600,'pt');
tiledlayout(2,4,'TileSpacing','tight','Padding','tight')

% Illlustrate neural and bold responses with different tau params
stimDeltaT = 0.1;
stimTimeTrial = (-3:stimDeltaT:30-stimDeltaT)';
stimProfile = [zeros(1,3/stimDeltaT),ones(1,12/stimDeltaT),zeros(1,18/stimDeltaT)]';
tauVals = [20,40,80,160];
hrf = returnFlobsVectors(stimDeltaT)*[0.86, 0.09, 0.01]';
hrf = hrf/sum(abs(hrf));


for tt = 1:length(tauVals)
    nexttile(tt)
    thisSignal = exp(-1/tauVals(tt)*stimTimeTrial).*stimProfile;
    thisSignal =  thisSignal / mean(thisSignal(thisSignal~=0));
    fit = conv2run(thisSignal,hrf,ones(size(stimProfile)));
    plot(stimTimeTrial,stimProfile,'-k');
    hold on
    plot(stimTimeTrial,thisSignal,'-r');
    plot(stimTimeTrial,fit,'-g');
    ylim([-0.1 1.5])
    a=gca();
    a.YTick = [0,1];
    title(sprintf('tau = %d',tauVals(tt)));
    axis off
end

nexttile([1,4])

for whichSet = 1:length(stimClassSetLabels)

    % Loop through subjects
    for ss = 1:length(subjectNames)

        % Load the results file for this subject
        fileName = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
        load(fileName,'results')

        % Load the results file for this subject
        fileName = fullfile(savePath,[subjectNames{ss} '_AvgV1_adaptMtSinai.mat']);
        load(fileName,'adaptResults')

        % Get the adaptation params
        pp(ss,:) = adaptResults.params(stimClassParamIdx{whichSet});

        x = log10(studiedFreqs) + (ss-2)*dirPlotShift;

        plot(x,pp(ss,:),subMarkers{ss},...
            'MarkerEdgeColor','none','MarkerFaceColor',plotColor{whichSet},...
            'MarkerSize',10)
        hold on

        % Report the difference in R2 value
        fprintf([subjectNames{ss} ' - R2 %2.3f, R2 adapt %2.3f\n'],max(results.R2),max(adaptResults.R2));

    end

    % Add median (across subject) values
    y = median(pp);
    yiqr = iqr(pp);
    for ff = 1:length(studiedFreqs)
        plot([log10(studiedFreqs(ff)) log10(studiedFreqs(ff))],[y(ff)-yiqr(ff)/2,y(ff)+yiqr(ff)/2],'-','Color',plotColor{whichSet},'LineWidth',4);
    end
    pH(whichSet) = plot(log10(studiedFreqs),y,'-','Color',lineColor{whichSet},'LineWidth',1.5);
    plot(log10(studiedFreqs),y,'o',...
        'MarkerFaceColor','none','MarkerEdgeColor',lineColor{whichSet},...
        'MarkerSize',15,'LineWidth',1.5)

end

% Clean up
legend(pH,stimClassSetLabels,'Location','northwest');
xlabel('Frequency [Hz]')
a=gca;
a.XTick = log10([2,4,8,16,32,64]);
a.XTickLabel = {'2','4','8','16','32','64'};
a.XTickLabelRotation = 0;
a.XMinorTick = 'off';
ylim([0 200]);
a.YTick = [0 50 100 150 200];
ylabel({'Time constant of','within-trial adaptation [secs]'})
xlabel('Stimulus frequency [Hz]')
box off

% Save the figure
%set(figHandle,'color','none');
fileName = fullfile(savePath,'WithinTrialAdaptation.pdf');
saveas(figHandle,fileName);

