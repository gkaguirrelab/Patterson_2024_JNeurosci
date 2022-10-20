% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','mriTemporalModel.mat');
load(loadPath,'mriTemporalModel');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriTemporalModel.meta.studiedFreqs;
studiedEccentricites = mriTemporalModel.meta.studiedEccentricites;
subjects = mriTemporalModel.meta.subjects;
stimulusDirections = mriTemporalModel.meta.stimulusDirections;
plotColor = mriTemporalModel.meta.plotColor;
nFixedParams = mriTemporalModel.meta.nFixedParams;
nFloatByEccParams = mriTemporalModel.meta.nFloatByEccParams;
nUniqueParams = mriTemporalModel.meta.nUniqueParams;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
cellClassOrder = {'midget','parasol','bistratified'};

for whichSub = 1:length(subjects)

    pMRI = mriTemporalModel.(subjects{whichSub}).pMRI;
    v1Y = mriTemporalModel.(subjects{whichSub}).data.v1Y;
    lgnY = mriTemporalModel.(subjects{whichSub}).data.lgnY;

    v1YFit = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams);
    lgnYFit = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,cellClassOrder,rgcTemporalModel,nUniqueParams,nFixedParams);

    for whichStim = 1:length(stimulusDirections)
        figure

        % Loop over eccentricities
        for ee=1:nEccs

            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;
            v1FitIndices = 1+(whichStim-1)*(nEccs*nFreqsForPlotting)+(ee-1)*(nFreqsForPlotting): ...
                (whichStim-1)*(nEccs*nFreqsForPlotting)+(ee-1)*(nFreqsForPlotting)+nFreqsForPlotting;

            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,v1Y(v1DataIndices),['o' plotColor{whichStim}]);
            hold on
            semilogx(freqsForPlotting,v1YFit(v1FitIndices),['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(studiedEccentricites(ee),2) '°']);
            ylim([-1 7]);
        end

        % Add the LGN response
        lgnDataIndices = 1+(whichStim-1)*(nFreqs): ...
            (whichStim-1)*(nFreqs)+nFreqs;
        lgnFitIndices = 1+(whichStim-1)*(nFreqsForPlotting): ...
            (whichStim-1)*(nFreqsForPlotting)+nFreqsForPlotting;

        subplot(2,4,8)
        semilogx(studiedFreqs,lgnY(lgnDataIndices),['o' plotColor{whichStim}]);
        hold on
        semilogx(freqsForPlotting,lgnYFit(lgnFitIndices),['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);
    end

    % Save the plot
    plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_MRIModelFit.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

    figure
    for whichCell = 1:length(cellClassOrder)
        % Plot the surround suppression index vs. eccentricity
        subplot(1,2,1)
        paramIndices = 1+nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams: ...
            nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams+nEccs;
        plot(log10(studiedEccentricites),pMRI(paramIndices),['*' plotColor{whichCell}]);
        hold on
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        subplot(1,2,2)
        semilogy(log10(studiedEccentricites),pMRI(paramIndices+nEccs),['*' plotColor{whichCell}]);
        hold on
        xlabel('Eccentricity [log deg]');
        ylabel('Gain parameter');
    end

        % Save the plot
    plotName = [subjects{whichSub} '_MRIModelParams.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

end
