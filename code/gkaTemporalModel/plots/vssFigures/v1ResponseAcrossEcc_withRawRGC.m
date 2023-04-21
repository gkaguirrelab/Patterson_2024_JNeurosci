% scriptCreatePlots

% Housekeeping
clear
%close all

% Properties of which model to plot
modelType = 'stimulus';
paramSearch = 'full';
freqsForPlotting = logspace(0,2,50);

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))),'data','temporalModelResults','v1',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

% Place to save figures
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_FullModel',modelType);

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
paramCounts = mriFullResultSet.meta.paramCounts;
cellClasses = {'midget','bistratified','parasol'};
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
nStims = 3;
myFreqs = logspace(log10(1),log10(100),101);
nFreqsForPlotting = length(myFreqs);
nCells = length(cellClasses);
subjectLineSpec = {'-','-'};

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.85 0.85 0.85],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = 0.4; % Transparency of the shaded error region

% Loop over subjects
for whichSub = 1:length(subjects)

    % Prepare the figures
    figHandles = figure('Renderer','painters');
    figuresize(400,800,'pt');
    tiledlayout(nEccs,3,'TileSpacing','tight','Padding','tight')

    % Get the model params and data
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);

    % Create a params vector that returns the unmodified RGC model output
    a = 0.2; % The gain on the RGC output
    pMRI = [1 120 1 0 a 0 a 0 a 1 0 0 0 0 0 0 a a a a a a 1 0 0 0 0 0 0 a a a a a a 1 0 0 0 0 0 0 a a a a a a];

    % Get the temporal model fits
    [~,v1YFitMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,paramCounts,modelType);

    % Loop over stimuli and plot
    for whichStim = 1:nStims

        % Get the average response across eccentricity
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Assemble the average V1 response for plotting later
            v1YThisEcc = v1Y(v1DataIndices);
            v1YThisEcclow = (v1Y(v1DataIndices)-v1YSEM(v1DataIndices));
            v1YThisEcchigh = (v1Y(v1DataIndices)+v1YSEM(v1DataIndices));

            % Now plot
            nexttile((ee-1)*nStims+stimOrder(whichStim));

            % Add a patch for the error
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ v1YThisEcclow, fliplr(v1YThisEcchigh) ],...
                plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha);
            hold on

            % Add the data symbols
            plot(log10(studiedFreqs),v1YThisEcc,...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',5,'MarkerEdgeColor','w','LineWidth',1);

            % Add the model fit
            plot(log10(freqsForPlotting),squeeze(v1YFitMatrix(whichStim,ee,:)),'-','Color',lineColor{stimOrder(whichStim)});

            % Add reference lines
            plot(log10([1 1]),[0 2],'-k');
            plot(log10([1 64]),[0 0],':k');

        end

    end

    % Clean up
    for ss=1:nStims
        for ee = 1:nEccs
            nexttile((ee-1)*nStims+ss);
            xlim(log10([0.5 150]))
            ylim([-1 7])
            a=gca;
            a.YTick = [0,2,4,6];
            a.YTickLabel = {'0','2','4','6'};
            a.XTick = log10([2,4,8,16,32,64]);
            a.XTickLabel = {'2','4','8','16','32','64'};
            a.XTickLabelRotation = 0;
            a.XMinorTick = 'off';
            a.YAxis.Visible = 'off';
            box off
            if ss>1
                a.XAxis.Visible = 'off';
                a.YAxis.Visible = 'off';
            end
        end
    end

    % Save the plots
    plotNamesPDF = [subjects{whichSub} '_v1ResponseAcrossEcc_withRawRGCModel.pdf' ];
    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end
