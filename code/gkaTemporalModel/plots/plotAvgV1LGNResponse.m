% scriptCreatePlots

% Housekeeping
clear
close all

modelType = 'stimulus';
paramSearches = {'full'};
paramLineStyle = {'-'};

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Where will we save the plots
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Create a figure
figure
figuresize(800,800,'pt');

% The order of plotting the stimulus directions
stimPlotOrder = [2 3 1];

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

    % Loop over the subjects
    for whichSub = 1:length(subjects)

        % Get some stuff from the result set
        pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
        pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
        v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
        v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
        v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
        lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
        lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
        lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

        % Obtain the fits
        [~,v1YFitMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);
        [~,lgnYFitMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);

        % Loop over stimulus directions
        for whichStim=1:length(stimulusDirections)

            % Set a subplot
            subIdx = (stimPlotOrder(whichStim)-1)*4+(whichSub-1)*2+2;
            subplot(3,4,subIdx);

            % Initialize variables to hold the average V1 response
            v1AvgY = zeros(1,length(studiedFreqs));
            v1YAvglow = zeros(1,length(studiedFreqs));
            v1YAvghigh = zeros(1,length(studiedFreqs));
            v1AvgFit = zeros(1,length(freqsForPlotting));

            % Loop over ecentricities
            for ee=1:length(studiedEccentricites)
                v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                    (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

                v1AvgY = v1AvgY + v1Y(v1DataIndices)./length(studiedEccentricites);
                v1YAvglow = v1YAvglow + v1Ylow(v1DataIndices)./length(studiedEccentricites);
                v1YAvghigh = v1YAvghigh + v1Yhigh(v1DataIndices)./length(studiedEccentricites);
                v1AvgFit = v1AvgFit + squeeze(v1YFitMatrix(whichStim,ee,:))/length(studiedEccentricites);
            end

            % Plot the data
                semilogx(studiedFreqs,v1AvgY,['.' plotColor{whichStim}],'MarkerSize',10);
                hold on

                % Add error bars
                X = [studiedFreqs fliplr(studiedFreqs)];
                Y = [v1YAvglow, fliplr(v1YAvghigh)];
                p = patch(X,Y,plotColor{whichStim});
                set(p,'edgecolor','none','facealpha',0.1);

            % Add the model fit
            semilogx(freqsForPlotting,v1AvgFit,[paramLineStyle{pp} plotColor{whichStim}]);

            % Clean up
            plot([16 16],[0 5.5],'--k');
            plot([1 128],[0 0],'-k');
            title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', V1']);
            ylim([-0.5 5.5]);
            xlim([1 128]);
            box off
            a=gca; a.XTick = [studiedFreqs(1) studiedFreqs(end)];
            a.XMinorTick = 'off';
            a.XTickLabel = arrayfun(@num2str, a.XTick, 'UniformOutput', 0);
            a.XTickLabelRotation = 0;
            a.YAxis.Visible = 'off';
            a.Color = 'none';

            if subIdx < 9
                a.XAxis.Visible = 'off';
            end


            % Now do the LGN
            subIdx = (stimPlotOrder(whichStim)-1)*4+(whichSub-1)*2+1;
            subplot(3,4,subIdx)

            % Add the LGN response
            lgnDataIndices = 1+(whichStim-1)*(nFreqs): ...
                (whichStim-1)*(nFreqs)+nFreqs;

                semilogx(studiedFreqs,lgnY(lgnDataIndices),['.' plotColor{whichStim}],'MarkerSize',10);
                hold on

                % Add error bars
                X = [studiedFreqs fliplr(studiedFreqs)];
                Y = [lgnYlow(lgnDataIndices), fliplr(lgnYhigh(lgnDataIndices))];
                p = patch(X,Y,plotColor{whichStim});
                set(p,'edgecolor','none','facealpha',0.1);

            % Model fit
            semilogx(freqsForPlotting,lgnYFitMatrix(whichStim,:),['-' plotColor{whichStim}]);

            % Clean up
            plot([16 16],[0 5.5],'--k');
            plot([1 128],[0 0],'-k');
            title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', LGN']);
            ylim([-0.5 5.5]);
            xlim([1 128]);
                        box off
            a=gca; a.XTick = [studiedFreqs(1) studiedFreqs(end)];
            a.XMinorTick = 'off';
            a.XTickLabel = arrayfun(@num2str, a.XTick, 'UniformOutput', 0);
            a.XTickLabelRotation = 0;
            a.YTick = 0:5;
            a.Color = 'none';

            if ~any([1 5 9]==subIdx) 
                        a.YAxis.Visible = 'off';
            end

            if subIdx < 9
                a.XAxis.Visible = 'off';
            end


        end % stims
    end % subjects
end % Loop over param searchs


% Save the plot
plotName = 'avgV1LGNRespose.pdf';
saveas(gcf,fullfile(savePath,plotName));