% scriptCreatePlots

% Housekeeping
clear
%close all

modelType = 'stimulus';
paramSearch = 'full';

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','v1',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

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
myFreqs = logspace(log10(1),log10(100),101);
nFreqsForPlotting = length(myFreqs);
nCells = length(cellClasses);
subjectLineSpec = {'-','-'};


% Params that control the plot appearance
spacing = 1;
stimOrder = [2 3 1];

plotColor={[0.85 0.85 0.85],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = [0,1,1];


for whichSub = 1:length(subjects)

    % Prepare the figures
    figHandles = figure('Renderer','painters');
    figuresize(200,400,'pt');
    set(gcf, 'Color', 'None');
    tiledlayout(3,1,'TileSpacing','none','Padding','tight')

    % Get the model params and data
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);

    % Loop over stimuli and plot
    for whichStim = 1:length(stimulusDirections)

        % Initialize variables to hold the average V1 response
        v1AvgY = zeros(1,length(studiedFreqs));
        v1YAvglow = zeros(1,length(studiedFreqs));
        v1YAvghigh = zeros(1,length(studiedFreqs));
        v1AvgFit = zeros(1,length(myFreqs));

        % Get the average response across eccentricity
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Assemble the average V1 response for plotting later
            v1AvgY = v1AvgY + v1Y(v1DataIndices)./length(studiedEccentricites);
            v1YAvglow = v1YAvglow + (v1Y(v1DataIndices)-v1YSEM(v1DataIndices))./length(studiedEccentricites);
            v1YAvghigh = v1YAvghigh + (v1Y(v1DataIndices)+v1YSEM(v1DataIndices))./length(studiedEccentricites);
        end

        % Now plot
        nexttile(stimOrder(whichStim));

        % Add a patch for the error
        

        % Add the error bars
        for bb=1:length(v1YAvglow)
            semilogx([studiedFreqs(bb) studiedFreqs(bb)], ...
                [v1YAvglow(bb), v1YAvghigh(bb)],...
                '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',2.0);
            hold on
        end

        % Add the data symbols
        semilogx(studiedFreqs,v1AvgY,...
            'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
            'MarkerSize',10,'MarkerEdgeColor','w','LineWidth',1);

        % Add reference lines
        semilogx([1 1],[0 2],'-k');
        semilogx([1 64],[0 0],':k');

    end

    % Clean up
    for ss=1:3
        nexttile(ss);
        xlim([0.5 150])
        ylim([-1 7])
        a=gca;
        a.YTick = [0,2,4,6];
        a.YTickLabel = {'0','2','4','6'};
        a.XTick = [2,4,8,16,32,64];
        a.XTickLabel = {'2','4','8','16','32','64'};
        a.XTickLabelRotation = 0;
        a.XMinorTick = 'off';
        a.YAxis.Visible = 'off';
        box off
        if ss<3
            a.XAxis.Visible = 'off';
            a.YAxis.Visible = 'off';
        end
    end

    % Save the plots
    plotNamesPDF = [subjects{whichSub} '_StackedV1Data.pdf' ];

    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end
