% scriptCreatePlots

% Housekeeping
clear
%close all

modelType = 'stimulus';
paramSearch = 'full';

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

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
subjectLineSpec = {'-','-'};


% Params that control the plot appearance
spacing = 3.5;
stimOrder = [2 3 1];

for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
        v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
        v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
        lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
        lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

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
            Y = [v1Y(v1DataIndices)-v1YSEM(v1DataIndices)+offset, fliplr(v1Y(v1DataIndices)+v1YSEM(v1DataIndices))+offset];
            p = patch(X,Y,plotColor{whichStim});
            set(p,'edgecolor','none','facealpha',0.2);

            % Add the model fit
            semilogx(fitFrequencies,squeeze(v1YFitMatrix(whichStim,ee,:))+offset,[subjectLineSpec{whichSub} plotColor{whichStim}],'LineWidth',2);

            % Add a refline
            semilogx([1 1],[offset offset+2],'-k');
            semilogx([1 100],[offset offset],':','Color',[0.5 0.5 0.5]);

            % Add a text label for the eccentricitiy
                text(80,offset+spacing/2,sprintf('%2.0f°',studiedEccentricites(ee)));

        end

        % Clean up
        semilogx([16 16],[0 (nEccs+0.5)*spacing],'--k')
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
