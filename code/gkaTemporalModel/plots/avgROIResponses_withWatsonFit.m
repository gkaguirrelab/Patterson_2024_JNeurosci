% scriptCreatePlots

% Housekeeping
clear
close all

% Properties of which model to plot
freqsForPlotting = logspace(0,2,50);

% Load the MRI data
mriData = loadMRIResponseData();

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Fixed features of the model
nCells = 3; nParams = 3;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];
interpFreqs = logspace(log10(1),log10(100),501);


% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [7 10 13]; % shifts each ttf down so they can be presented tightly on the same figure

roiNames = {'lgn','v1_avg','v2v3_avg'};
nROIs = length(roiNames);
nAcqs = 12;

% fmincon Options. Indicate that the objective function is deterministic,
% and handle verbosity
options = optimoptions('fmincon');
options.Display = 'none';

optionsBADS.UncertaintyHandling = 0;
optionsBADS.Display = 'off';

% Set some bounds
LB = [0 1 0.5 0.5];
UB = [5 10 3 3];
p0 = [1.5 5 1 1.5];

% Loop over subjects
for whichSub = 1:length(subjects)

    % Prepare the figure
    figHandles = figure('Renderer','painters');
    figuresize(600,400,'pt');
    tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

    % Loop over ROIS
    for rr = 1:nROIs

        % Load the data
        Y = zeros(nStims,length(studiedFreqs));
        Ysem = zeros(nStims,length(studiedFreqs));
        for whichStim = 1:nStims
            thisMatrix = mriData.(subjects{whichSub}).(stimulusDirections{whichStim}).(roiNames{rr});
            Y(whichStim,:) = mean(thisMatrix);
            Ysem(whichStim,:) = std(thisMatrix)/sqrt(nAcqs);
        end

        % Loop over stimuli and plot
        for whichStim = 1:nStims

            % Assemble the data
            YThisROI = squeeze(Y(whichStim,:));
            YsemThisROI = squeeze(Ysem(whichStim,:));
            YThisROIlow = YThisROI - YsemThisROI;
            YThisROIhigh = YThisROI + YsemThisROI;

            % The weighted objective
            myObj = @(p) norm( (1./YsemThisROI) .* ( YThisROI - watsonTemporalModel(p,studiedFreqs)));

            % Fit it
            p = fmincon(myObj,p0,[],[],[],[],LB,UB,[],options);

            % Get the fitted response
            yFit = watsonTemporalModel(p,freqsForPlotting);

            % Select the plot of the correct stimulus direction
            nexttile(stimOrder(whichStim));

            % Add a patch for the error
            patch(...
                [log10(studiedFreqs),fliplr(log10(studiedFreqs))],...
                [ YThisROIlow-shift_ttf(rr), fliplr(YThisROIhigh)-shift_ttf(rr) ],...
                plotColor{whichStim},'EdgeColor','none','FaceColor',plotColor{stimOrder(whichStim)},'FaceAlpha',faceAlpha);
            hold on

            % Add the data symbols, using reversed markers for values below
            % zero
            idx = YThisROI > 0;
            plot(log10(studiedFreqs(idx)),YThisROI(idx)-shift_ttf(rr),...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);
            idx = YThisROI < 0;
            plot(log10(studiedFreqs(idx)),YThisROI(idx)-shift_ttf(rr),...
                'o','MarkerFaceColor','w',...
                'MarkerSize',6,'MarkerEdgeColor',lineColor{stimOrder(whichStim)},'LineWidth',1);

            % Add the model fit
            plot(log10(freqsForPlotting),yFit-shift_ttf(rr),...
                ['-' lineColor{stimOrder(whichStim)}],...
                'LineWidth',2);

            % Add reference lines
            if rr==1 && whichStim == 3
                plot(log10([1 1]),[0 2],'-k');
            end
            plot(log10([1 2]),[0 0]-shift_ttf(rr),':k');
            plot(log10([50 100]),[0 0]-shift_ttf(rr),':k');

        end

    end

    % Clean up
    for ss=1:nStims
        nexttile(ss);
        xlim(log10([0.5 150]))
        ylim([-14 5])
        xlabel('Frequency [Hz]')
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
    plotNamesPDF = [subjects{whichSub} '_avgROIResponses_withWatsonModel.pdf' ];
    saveas(figHandles,fullfile(savePath,plotNamesPDF));

end





