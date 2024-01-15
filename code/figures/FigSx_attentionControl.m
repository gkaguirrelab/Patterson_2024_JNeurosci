% This script creates an illustrative figure of the average time-series
% data and model fit over all of area V1.

clear
close all
rng(1);

%% Analysis properties
nBoots = 250;

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names and stimulus directions
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
shortNames = {'gka','asb','cgp'};
stimulusDirections = {'LminusM','S','LMS'};
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};

nDirections = length(stimulusDirections);

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
stimOrder = [2 3 1];
lineStyles ={'-','--'};
interpFreqs = logspace(log10(1),log10(100),501);

figHandle = figure('Renderer','painters');
figuresize(600,300,'pt');

% Loop through the subjects
for ss = 1:length(subjectNames)

    subplot(1,3,ss)

    for aa = 1:2

        % Load the results file for this subject
        switch aa
            case 1
                filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
            case 2
                filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_attentionControl_results.mat']);
        end

        load(filePath,'results')

        % This analysis was conducted upon the average V1 signal. We just need
        % a single index for a vertex that was within V1.
        goodIdx = find(~isnan(results.fVal),1);

        % Grab the stimLabels
        stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

        % Loop over the stimuli
        for whichStim = 3:3

            % Loop through the stimuli and frequencies and obtain the
            % data for the trials that contained attention events and
            % those that did not
            rawVals = [];
            adjustedVals = [];

            for ff = 1:length(allFreqs)
                % By searching for stimulus labels that begin with "f", we
                % exclude those trials that contained an attention event
                % under case 2 (when we have loaded the "attention control"
                % result set.
                subString = sprintf(['f%dHz_' stimulusDirections{whichStim}],allFreqs(ff));
                idx = find(startsWith(stimLabels,subString));
                params = results.params(goodIdx,idx);
                % Check here if we were unlucky enough to not have a trial
                % type that occurred without an attention event
                if any(params(:)==0)
                    foo=1;
                end
                rawVals(ff,:) = params;
            end

            % Adjust the values for the zero frequency
            for ff = 2:length(allFreqs)
                adjustedVals(ff-1,:) = rawVals(ff,:)-rawVals(1,:);
            end

            % Get the median across (bootstrap resampled) acquisitions
            for bb = 1:nBoots
                bootIdx = datasample(1:length(idx),length(idx));
                Yboot(bb,:) = mean(adjustedVals(:,bootIdx),2)';
            end

            Y = median(Yboot)';
            YIQR = iqr(Yboot)';

            % Get the Watson fit
            [p,fVal,~,yFitInterp] = fitWatsonModel(Y',1./YIQR',studiedFreqs);

            % Plot
            lowY = (Y - YIQR/2)';
            highY = (Y + YIQR/2)';

            % Add a patch for the error
            if aa == 1
                patch(...
                    [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                    [ lowY, fliplr(highY) ],...
                    plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha/aa);

                hold on

                % Add the data symbols, using reversed markers for values below
                % zero
                idx = Y > 0;
                plot(log10(studiedFreqs(idx)),Y(idx),...
                    'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                    'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);
                idx = Y < 0;
                plot(log10(studiedFreqs(idx)),Y(idx),...
                    'o','MarkerFaceColor','w',...
                    'MarkerSize',6,'MarkerEdgeColor',lineColor{stimOrder(whichStim)},'LineWidth',1);

            end

            % Add the model fit
            ph(aa) = plot(log10(interpFreqs),yFitInterp,...
                [lineStyles{aa} lineColor{stimOrder(whichStim)}],...
                'LineWidth',2);

        end % with and without attention event trials

    end % stimuli

    title(subjectNames{ss})
    ylim([0 5]);
    a=gca;
    a.YTick = [0,5];
    a.YTickLabel = {'0','5'};
    a.XTick = log10([2,4,8,16,32,64]);
    a.XTickLabel = {'2','4','8','16','32','64'};
    a.XTickLabelRotation = 0;
    a.XMinorTick = 'off';


    if ss==1
        legend(ph,{'all trials','exclude attention'})
    end

end


% Save the figure
set(figHandle,'color','none');
fileName = fullfile(savePath,'LMS_V1_TTF_withAndwithoutAttentionTrials.pdf');
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
