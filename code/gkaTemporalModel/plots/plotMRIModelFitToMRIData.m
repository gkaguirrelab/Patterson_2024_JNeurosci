% scriptCreatePlots

% Housekeeping
clear
close all

% Model and search type
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

% Where will we save the plots
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
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
nCells = length(cellClasses);
subjectLineSpec = {'-','--'};

relGainFig = figure();

% Loop over subjects and make plots
for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

    % Calculate the gain ratios
    gainVals = {};
    for whichStim = 1:length(stimulusDirections)
        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs + 1;
        gainVals(whichStim) = {mriFullResultSet.(subjects{whichSub}).pMRI(:,startIdx:startIdx+nEccs-1)};
    end

    % Plot the gain ratios    
    figure(relGainFig)
    subplot(1,2,2)
    semilogy(log10(studiedEccentricites),mean(gainVals{1}./gainVals{3}),[subjectLineSpec{whichSub} 'r']);
    hold on
    semilogy(log10(studiedEccentricites),mean(gainVals{2}./gainVals{3}),[subjectLineSpec{whichSub} 'b']);
    semilogy([-0.5 2],[1 1],':k');
    xlim([-0.5 2]);
    ylim([10^-1 10^2]);
    xlabel('Eccentricity [log deg]');
    ylabel('Relative gain parameter');
    box off


    % Get the temporal model fits
    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);
    [~,lgnYFitMatrix,lgnRFMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);

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
            Y = [v1Y(v1DataIndices)-v1YSEM(v1DataIndices), fliplr(v1Y(v1DataIndices)+v1YSEM(v1DataIndices))];
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
        Y = [lgnY(lgnDataIndices)-lgnYSEM(lgnDataIndices), fliplr(lgnY(lgnDataIndices)+lgnYSEM(lgnDataIndices))];
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

        % Show the synapse LP filter
        subplot(2,4,4)
        synapseSecondOrderFc = pMRI(2);
        synapseSecondOrderQ = 0.45;

        rf = stageSecondOrderLP(synapseSecondOrderFc,synapseSecondOrderQ);
        myFreqs = logspace(log10(0.5),log10(100),101);
        ttfComplex = double(subs(rf,myFreqs));
        gainVals = abs(ttfComplex);
        semilogx(myFreqs,gainVals,'-k');

        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;

        % Save the plot
        plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_' paramSearch '_MRIModelFit.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

    end

    % Get the gain for achromatic stimuli at the fovea
    whichStim = 3;
    referenceIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs +1;
    gainReference = pMRI(referenceIdx);

    % Plot the params across eccentricity
    figure
    for whichStim = 1:length(stimulusDirections)

        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed;
        v1SurroundIndex = pMRI(startIdx+1:startIdx+nEccs);
        v1SurroundIndexSEM = pMRISEM(startIdx+1:startIdx+nEccs);
        v1Gain = pMRI(startIdx+1+nEccs:startIdx+nEccs+nEccs);
        v1GainSEM = pMRISEM(startIdx+1+nEccs:startIdx+nEccs+nEccs);

        % Plot the surround suppression index vs. eccentricity
        subplot(1,2,1)
        plot(log10(studiedEccentricites),v1SurroundIndex,['o' plotColor{whichStim}]);
        hold on
        plot(log10(studiedEccentricites),v1SurroundIndex,[subjectLineSpec{whichSub} plotColor{whichStim}]);

        % Add error bars
        X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
        Y = [v1SurroundIndex-v1SurroundIndexSEM, fliplr(v1SurroundIndex+v1SurroundIndexSEM)];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.1);

        % Add the lgn
        idx = paramCounts.unique+(whichStim-1)*paramCounts.lgn+1;
        semilogy(0,pMRI(idx),['*' plotColor{whichStim}]);
        semilogy([0 0],[pMRI(idx)-pMRISEM(idx), pMRI(idx)+pMRISEM(idx)],['-' plotColor{whichStim}]);

        % Clean up
        semilogy([-0.5 2],[0 0],':k');
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        xlim([-0.5 2]);
        ylim([-1 1]);

        % Plot the gain vs. eccentricity
        subplot(1,2,2)
        semilogy(log10(studiedEccentricites),v1Gain/gainReference,['o' plotColor{whichStim}]);
        hold on
        semilogy(log10(studiedEccentricites),v1Gain/gainReference,[subjectLineSpec{whichSub} plotColor{whichStim}]);

        % Add error bars
        Y = [v1Gain/gainReference-v1GainSEM/gainReference, fliplr(v1Gain/gainReference+v1GainSEM/gainReference)];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.1);

        % Add the lgn and lgn error bars
        idx = paramCounts.unique+(whichStim-1)*paramCounts.lgn+2;
        semilogy(0,pMRI(idx)/gainReference,['*' plotColor{whichStim}]);
        semilogy([0 0],[pMRI(idx)/gainReference-pMRISEM(idx)/gainReference, pMRI(idx)/gainReference+pMRISEM(idx)/gainReference],['-' plotColor{whichStim}]);

        % Clean up
        semilogy([-0.5 2],[1 1],':k');
        xlim([-0.5 2]);
        ylim([10^-1 10^2]);
        xlabel('Eccentricity [log deg]');
        ylabel('Relative gain parameter');

    end

    % Save the parameter plot
    plotName = [subjects{whichSub} '_' paramSearch '_MRIModelParams.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));


    %
    %     % For one subject and one eccentricity, plot the lgn and V1 IRFs
    %     ee=4;
    %     if whichSub == 1
    %         figHandleLGN = figure('Position',  [100, 100, 200, 400]);
    %         figHandleV1 = figure('Position',  [100, 100, 200, 400]);
    %         for ss=1:length(stimulusDirections)
    %             plotRF(lgnRFMatrix(ss,ee),figHandleLGN,['-' plotColor{ss}]);
    %             plotRF(v1RFMatrix(ss,ee),figHandleV1,['-' plotColor{ss}]);
    %         end
    %
    %         % Save the plot
    %         plotName = [subjects{whichSub} '_' paramSearch '_lgnIRFs.pdf' ];
    %         saveas(figHandleLGN,fullfile(savePath,plotName));
    %         plotName = [subjects{whichSub} '_' paramSearch '_v1IRFs.pdf' ];
    %         saveas(figHandleV1,fullfile(savePath,plotName));
    %     end

end

    % Save across-subject parameter plot
    plotName = [paramSearch '_RelativeGainAcrossSubjects.pdf' ];
    saveas(relGainFig,fullfile(savePath,plotName));

