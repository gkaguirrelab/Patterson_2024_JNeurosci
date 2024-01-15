% This script creates an illustrative figure of the average time-series
% data and model fit over all of area V1.

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

    % Properties of the stimuli
    stimLabels = results.model.opts{6};

    % Plot the average time series and model fit
    for whichStim=1:3
        data = results.data.avgSignal(1+(whichStim-1)*672:672+(whichStim-1)*672);
        fit = results.data.avgModelFit(1+(whichStim-1)*672:672+(whichStim-1)*672);
        residual = data - fit;
        saveResid(ss,whichStim,:) = residual;
        saveFit(ss,whichStim,:) = fit;
    end
end

fit = squeeze(mean(mean(saveFit)));
residual = squeeze(mean(mean(saveResid)));

figure

plot(fit,'-r');
hold on
plot(residual,'-k');
