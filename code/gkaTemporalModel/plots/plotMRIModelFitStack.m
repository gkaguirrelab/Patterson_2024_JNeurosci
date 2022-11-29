% scriptCreatePlots

% Housekeeping
clear
%close all

modelType = 'stimulus';

zeroSurroundEffect = false;

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
load(fullfile(loadPath,'mriFullResultSet.mat'),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_FullModel',modelType);

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
paramCounts = mriFullResultSet.meta.paramCounts;
cellClasses = {'midget','bistratified','parasol'};
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
fitFrequencies = logspace(0,2,50);
nFreqsForPlotting = length(fitFrequencies);
nCells = length(cellClasses);


% Params that control the plot appearance
spacing = 3.5;
stimOrder = [2 3 1];

for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
    lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
    lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
    lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

    if zeroSurroundEffect
%        pMRI([16:21, 29:34, 42:47]) = 0;
    end

    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,fitFrequencies,rgcTemporalModel,paramCounts,modelType);

    figure
    figuresize(500,300,'pt');

    for whichStim = 1:length(stimulusDirections)

        subplot(1,3,stimOrder(whichStim));

        % Loop over eccentricities
        for ee=1:nEccs

            offset = 1+(nEccs*spacing) - (ee)*spacing;

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Show the data itself
            semilogx(studiedFreqs,v1Y(v1DataIndices)+offset,...
                '.','Color',plotColor{whichStim},'MarkerSize',10);
            hold on

            % Add error bars
            X = [studiedFreqs fliplr(studiedFreqs)];
            Y = [v1Ylow(v1DataIndices), fliplr(v1Yhigh(v1DataIndices))]+offset;
            p = patch(X,Y,plotColor{whichStim});
            set(p,'edgecolor','none','facealpha',0.1);

            % Add the model fit
            semilogx(fitFrequencies,squeeze(v1YFitMatrix(whichStim,ee,:))+offset,['-' plotColor{whichStim}]);

            % Add a refline
            semilogx([1 1],[offset offset+2],'-k');
            semilogx([1 100],[offset offset],':','Color',[0.5 0.5 0.5]);

        end

        % Clean up
        title([stimulusDirections{whichStim} ', ' subjects{whichSub} ]);
        xlim([0.5 100])
        ylim([0 23])
        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;
        a.XMinorTick = 'off';
        a.YAxis.Visible = 'off';
        box off
    end

    % Save the plot
        plotName = [subjects{whichSub} '_StackedModelFits.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

end
