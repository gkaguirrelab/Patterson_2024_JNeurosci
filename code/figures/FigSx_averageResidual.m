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
localSaveDir = getpref('Patterson_2024_JNeurosci','localSaveDir');

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
    figure
    resid = results.data.datats-results.data.modelts;
    stimMat = cell2mat(results.model.inputs{2});
stimLabels = results.model.opts{4};

    freqSet = {'2Hz','4Hz','8Hz','16Hz','32Hz','64Hz'};

    hrfDelay = 5;

    for ff=1:length(freqSet)
        idx=find(contains(stimLabels,freqSet{ff}));

        trialResp = [];
        for ii=1:length(idx)
            vec = stimMat(idx(ii),:);
            starts = find(diff(vec)>0.99);
            for ss=1:length(starts)

                startIdx = starts(ss)+hrfDelay;
                endIdx = starts(ss)+14+hrfDelay;
                if endIdx<=length(resid)
                trialResp(end+1,:) = resid(startIdx:endIdx);
                end
            end
        end
        plot(mean(trialResp));
        hold on
    end
    ylim([-1 1]);

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
