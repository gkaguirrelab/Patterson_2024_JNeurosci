% scriptCreatePlots

% Housekeeping
clear
%close all

modelType = 'stimulus';
paramSearch = 'full';

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))),'data','temporalModelResults','v1',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

% Place to save figures
savePath = '/Users/carlynpattersongentile/Aguirre-Brainard Lab Dropbox/Carlyn Patterson Gentile/Patterson_2021_EccentricityFlicker/VSS 2023/figures/';

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
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [0 3 6 9 11 13]; % shifts each ttf down so they can be presented tightly on the same figure

% Loop over subjects
for whichSub = 1:length(subjects)

    % Prepare the figures
    figHandles = figure('Renderer','painters');
    figuresize(600,400,'pt');
    tiledlayout(1,3,'TileSpacing','tight','Padding','tight')
    

    % Get the model params and data
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);

    % Loop over stimuli and plot
    for whichStim = 1:nStims

        % Get the average response across eccentricity
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Assemble the average V1 response for plotting later
            v1YThisEcc = v1Y(v1DataIndices);
            v1YThisEcclow = (v1Y(v1DataIndices)-v1YSEM(v1DataIndices)-shift_ttf(ee));
            v1YThisEcchigh = (v1Y(v1DataIndices)+v1YSEM(v1DataIndices)-shift_ttf(ee));

            % select the correct stimulus plot
            nexttile(stimOrder(whichStim));

            % Add a patch for the error
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ v1YThisEcclow, fliplr(v1YThisEcchigh) ],...
                plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha);
            hold on

            % Add the data symbols
            plot(log10(studiedFreqs),v1YThisEcc-shift_ttf(ee),...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);

            % Add reference lines
            if ee==1 && whichStim == 3
                plot(log10([1 1]),[0 2],'-k');
            end
                plot(log10([1 2]),[0 0]-shift_ttf(ee),':k');

        end

    end

    % Clean up
    for ss=1:nStims
        nexttile(ss);
        xlim(log10([0.5 150]))
        ylim([-14 5])
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
        end
    end

    % Save the plots
    plotNamesPDF = [subjects{whichSub} '_v1ResponseAcrossEcc.pdf' ];
    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end