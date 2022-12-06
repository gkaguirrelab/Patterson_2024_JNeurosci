% scriptCreatePlots

% Housekeeping
clear
%close all

modelType = 'stimulus';
paramSearches = {'gainOnly','zeroSurroundIndex','full'};
paramLineStyle = {':','-.','-'};

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Where will we save the plots
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Create a figure
figure
figuresize(200,300,'pt');

whichSub = 2;
whichEcc = 4;
whichStim = 3;

% Loop over the param searches
for pp = 1:length(paramSearches)

    % Load the MRI temporal model
    loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults',modelType);
    load(fullfile(loadPath,['mriFullResultSet_' paramSearches{pp} '.mat']),'mriFullResultSet');

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
    freqsForPlotting = logspace(0,2,50);
    nFreqsForPlotting = length(freqsForPlotting);
    nCells = length(cellClasses);

    % Grabt the materials for this subject

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;

    [~,v1YFitMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);

    % Get the lgn fit with the separate parasol and midget
    % contributions to LMS
    if pp==1
        [~,v1YFitMatrixMidgetOnly] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType,{'midget'});
        [~,v1YFitMatrixParasolOnly] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType,{'parasol'});
    end


    % The indices of the data to be plotted in the big vector
    v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(whichEcc-1)*(nFreqs): ...
        (whichStim-1)*(nEccs*nFreqs)+(whichEcc-1)*(nEccs)+nFreqs;

    % Show the data itself
    if pp==1
        semilogx(studiedFreqs,v1Y(v1DataIndices),['o' plotColor{whichStim}]);
        hold on

        % Add error bars
        X = [studiedFreqs fliplr(studiedFreqs)];
        Y = [v1Ylow(v1DataIndices), fliplr(v1Yhigh(v1DataIndices))];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.1);

        if whichStim == 3
            semilogx(freqsForPlotting,squeeze(v1YFitMatrixMidgetOnly(whichStim,whichEcc,:)),':r');
            semilogx(freqsForPlotting,squeeze(v1YFitMatrixParasolOnly(whichStim,whichEcc,:)),':','Color',[0.5 0.5 0.5]);
        end
    end

    % Add the model fit
    semilogx(freqsForPlotting,squeeze(v1YFitMatrix(whichStim,whichEcc,:)),[paramLineStyle{pp} plotColor{whichStim}]);

    % Clean up
    refline(0,0);
    title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', V1']);
    ylim([-1 7]);
    a=gca; a.XTick = studiedFreqs;
    a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
    a.XTickLabelRotation = 0;

end % Loop over param searchs


% Save the plot
plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_modelStagesIllustration.pdf' ];
saveas(gcf,fullfile(savePath,plotName));
