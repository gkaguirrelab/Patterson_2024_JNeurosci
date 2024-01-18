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
subjects = {'gka','asb','cgp'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};

% Color map
cmap = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';

% Loop through the subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(savePath,[subjectNames{ss} '_avgV1_carryOverMtSinai.mat']);
    load(filePath,'carryOverResults')

    % reshape the carry-over beta values into a matrix
    paramStartIdx = find(startsWith(carryOverResults.model.opts{4},'co_f0->f0_'));
    for dd = 1:3
        coParams = carryOverResults.params(paramStartIdx(dd):paramStartIdx(dd)+48);
        coLabels = carryOverResults.model.opts{4}(paramStartIdx(dd):paramStartIdx(dd)+48);
        for ii=1:49
            [rr,cc] = ind2sub([7 7],ii);
            coMatrix{ss,dd}(rr,cc)=coParams(ii);
        end
    end

    % Plot the matrices
    figHandle = figure();
    figuresize(600,300,'pt');
    t = tiledlayout(1,3);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    for dd=1:3
        nexttile(dd);
        im = coMatrix{ss,dd};
        imagesc(im);
        colormap(cmap)
        axis square
        title([subjects{ss} ' - ' directions{dd} ]);
        drawnow

    end
end
