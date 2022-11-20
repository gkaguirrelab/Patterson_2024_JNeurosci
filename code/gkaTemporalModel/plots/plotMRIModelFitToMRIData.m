% scriptCreatePlots

% Housekeeping
clear
%close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults');
load(fullfile(loadPath,'mriFullResultSet.mat'),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_FullModel');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
paramCounts = mriFullResultSet.meta.paramCounts;
postReceptoralPaths = mriFullResultSet.meta.postReceptoralPaths;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
cellClasses = {'midget','bistratified','parasol'};
nCells = length(cellClasses);

for whichSub = 1:length(subjects)


    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
    lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
    lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
    lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts);
    [~,lgnYFitMatrix,lgnRFMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,freqsForPlotting,rgcTemporalModel,paramCounts);

    for whichStim = 1:length(stimulusDirections)
        figure
        figuresize(600,300,'pt');

        % Loop over eccentricities
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Create a subplot
            subplot(2,4,ee+(ee>3))

            % Show the data itself
            semilogx(studiedFreqs,v1Y(v1DataIndices),['o' plotColor{whichStim}]);
            hold on

            % Add error bars
            X = [studiedFreqs fliplr(studiedFreqs)];
            Y = [v1Ylow(v1DataIndices), fliplr(v1Yhigh(v1DataIndices))];
            p = patch(X,Y,plotColor{whichStim});
            set(p,'edgecolor','none','facealpha',0.1);

            % Add the model fit
            semilogx(freqsForPlotting,squeeze(v1YFitMatrix(whichStim,ee,:)),['-' plotColor{whichStim}]);

            % Clean up
            refline(0,0);
            title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(studiedEccentricites(ee),2) '°']);
            ylim([-1 7]);
            a=gca; a.XTick = studiedFreqs;
            a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
            a.XTickLabelRotation = 0;

        end

        % Add the LGN response
        lgnDataIndices = 1+(whichStim-1)*(nFreqs): ...
            (whichStim-1)*(nFreqs)+nFreqs;

        % Show the data
        subplot(2,4,8)
        semilogx(studiedFreqs,lgnY(lgnDataIndices),['o' plotColor{whichStim}]);
        hold on

        % patch error bars
        X = [studiedFreqs fliplr(studiedFreqs)];
        Y = [lgnYlow(lgnDataIndices), fliplr(lgnYhigh(lgnDataIndices))];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.1);

        % Model fit
        semilogx(freqsForPlotting,lgnYFitMatrix(whichStim,:),['-' plotColor{whichStim}]);

        % Clean up
        refline(0,0);
        title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-1 7]);
        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;

        % Show the two filters
        subplot(2,4,4)
        lgnSecondOrderFc = pMRI(2);
        lgnSecondOrderQ = pMRI(3);
        v1SecondOrderFc = pMRI(4);
        v1SecondOrderQ = pMRI(5);

        rf = stageSecondOrderLP(lgnSecondOrderFc,lgnSecondOrderQ);
        myFreqs = logspace(log10(0.5),log10(100),101);
        ttfComplex = double(subs(rf,myFreqs));
        gainVals = abs(ttfComplex);
        semilogx(myFreqs,gainVals,'-k');
        hold on
        rf = stageSecondOrderLP(v1SecondOrderFc,v1SecondOrderQ);
        myFreqs = logspace(log10(0.5),log10(100),101);
        ttfComplex = double(subs(rf,myFreqs));
        gainVals = abs(ttfComplex);
        semilogx(myFreqs,gainVals,'-m');

        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;

        % Save the plot
        plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_MRIModelFit.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

    end

    % Plot the params across eccentricity
    figure
    for whichStim = 1:length(stimulusDirections)

        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed;
        v1SurroundIndex = pMRI(startIdx+1:startIdx+nEccs);
        v1Gain = pMRI(startIdx+1+nEccs:startIdx+nEccs+nEccs);

        % Plot the surround suppression index vs. eccentricity
        subplot(1,2,1)
        plot(log10(studiedEccentricites),v1SurroundIndex,['o' plotColor{whichStim}]);
        hold on
        plot(log10(studiedEccentricites),v1SurroundIndex,['-' plotColor{whichStim}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Plot the gain index vs. eccentricity
        subplot(1,2,2)
        semilogy(log10(studiedEccentricites),v1Gain,['o' plotColor{whichStim}]);
        hold on
        semilogy(log10(studiedEccentricites),v1Gain,['-' plotColor{whichStim}]);
        xlim([-0.1 2]);
        xlabel('Eccentricity [log deg]');
        ylabel('Gain parameter');
    end

    % Save the plot
    plotName = [subjects{whichSub} '_MRIModelParams.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

    % For one subject and one eccentricity, plot the lgn and V1 IRFs
    ee=4;
    if whichSub == 1
            figHandleLGN = figure('Position',  [100, 100, 200, 400]);
    figHandleV1 = figure('Position',  [100, 100, 200, 400]);
        for ss=1:length(stimulusDirections)
            plotRF(lgnRFMatrix(ss,ee),figHandleLGN,['-' plotColor{ss}]);
            plotRF(v1RFMatrix(ss,ee),figHandleV1,['-' plotColor{ss}]);
        end

    % Save the plot
    plotName = [subjects{whichSub} '_lgnIRFs.pdf' ];
    saveas(figHandleLGN,fullfile(savePath,plotName));
    plotName = [subjects{whichSub} '_v1IRFs.pdf' ];
    saveas(figHandleV1,fullfile(savePath,plotName));
    end
    
end
