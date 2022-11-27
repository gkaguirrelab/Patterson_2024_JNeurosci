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
plotColor = {[1 0 0],[0 0 1],[1 1 1]};

%fitEccentricities = logspace(log10(0.7031),log10(90),15);
fitEccentricities = studiedEccentricites;

spacing = 3;

for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
    lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
    lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
    lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,fitEccentricities,fitFrequencies,rgcTemporalModel,paramCounts,modelType);

    figure
    figuresize(600,300,'pt');

    for whichStim = 1:length(stimulusDirections)

        subplot(1,3,whichStim);

        % Loop over eccentricities
        for ee=nEccs:-1:1

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Create a patch plot of the fit
            X = [fitFrequencies fliplr(fitFrequencies)];
            Y = [squeeze(v1YFitMatrix(whichStim,ee,:))'+ee*spacing, repmat(ee*spacing,1,length(fitFrequencies))];
            p = patch(X,Y,plotColor{whichStim}.*0.5);
            hold on
            set(p,'edgecolor','w','facealpha',1,'LineWidth',3);


            plot(studiedFreqs,v1Y(v1DataIndices)+ee*spacing,...
                '.','Color',plotColor{whichStim},'MarkerSize',10);
        end


        % Clean up
        set(gca,'XScale','log');
        title([stimulusDirections{whichStim} ', ' subjects{whichSub} ]);
        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;

    end

    % Save the plot
    %     plotName = [subjects{whichSub} '_MRIModelParams.pdf' ];
    %     saveas(gcf,fullfile(savePath,plotName));

end

%