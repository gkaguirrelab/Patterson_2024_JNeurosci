% scriptCreatePlots

% Housekeeping
clear
%close all

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
    figHeights = [600,200];
    viewAngles = [70,20];
    for ff=1:2
        figHandles{ff} = figure('Renderer','painters');
        figuresize(600,figHeights(ff),'pt');
        set(gcf, 'Color', 'None');
        tiledlayout(1,3,'TileSpacing','none','Padding','tight')
    end

    % Get the model params and data
    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

    % Get the model fits
    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,myFreqs,rgcTemporalModel,paramCounts,modelType);
    [~,lgnYFitMatrix] = assembleLGNResponse(pMRI,cellClasses,stimulusDirections,myFreqs,rgcTemporalModel,paramCounts,modelType);

    % Loop over stimuli and plot
    for whichStim = 1:length(stimulusDirections)

        figure(figHandles{1});

        % Initialize variables to hold the average V1 response
        v1AvgY = zeros(1,length(studiedFreqs));
        v1YAvglow = zeros(1,length(studiedFreqs));
        v1YAvghigh = zeros(1,length(studiedFreqs));
        v1AvgFit = zeros(1,length(myFreqs));

        % Loop over eccentricities
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Assemble the average V1 response for plotting later
            v1AvgY = v1AvgY + v1Y(v1DataIndices)./length(studiedEccentricites);
            v1YAvglow = v1YAvglow + (v1Y(v1DataIndices)-v1YSEM(v1DataIndices))./length(studiedEccentricites);
            v1YAvghigh = v1YAvghigh + (v1Y(v1DataIndices)+v1YSEM(v1DataIndices))./length(studiedEccentricites);
            v1AvgFit = v1AvgFit + squeeze(v1YFitMatrix(whichStim,ee,:))'/length(studiedEccentricites);

            % First plot the fit
            yVals = squeeze(v1YFitMatrix(whichStim,ee,:))';

            % Get the subplot and define the offset
            nexttile(stimOrder(whichStim));
            set(gca, 'XScale', 'log')

            % Create a patch for the response
            X = [myFreqs fliplr(myFreqs)];
            Y = repmat(ee*spacing,size(X));
            Z = [yVals, zeros(size(yVals))];
            p = fill3(X,Y,Z,plotColor{stimOrder(whichStim)});
            set(p,'edgecolor','none','facealpha',1);
            hold on

            % Add a line at the edge of the color patch
            plot3(myFreqs,repmat(ee*spacing,size(myFreqs)),yVals,...
                '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.5);

            % Add the error bars
            yVals = v1Y(v1DataIndices);
            ySEM = v1YSEM(v1DataIndices);
            for bb=1:length(ySEM)
                plot3([studiedFreqs(bb) studiedFreqs(bb)],[ee*spacing ee*spacing], ...
                    [yVals(bb)-ySEM(bb), yVals(bb)+ySEM(bb)],...
                    '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.0);
            end

            % Add the data symbols
            plot3(studiedFreqs,repmat(ee*spacing,size(studiedFreqs)),yVals,...
                'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
                'MarkerSize',5,'MarkerEdgeColor','w','LineWidth',0.5);

            % Add a text label for the eccentricitiy
            text(150,ee*spacing,0,sprintf('%2.0f°',studiedEccentricites(ee)),"HorizontalAlignment","left");

        end

        % Now plot average V1 and LGN responses
        figure(figHandles{2});
        nexttile(stimOrder(whichStim));

        % LGN
        % Create a patch for the response
        X = [myFreqs fliplr(myFreqs)];
        Y = repmat(1,size(X));
        Z = [lgnYFitMatrix(whichStim,:), zeros(size(myFreqs))];
        p = fill3(X,Y,Z,plotColor{stimOrder(whichStim)});
        set(p,'edgecolor','none','facealpha',1);
        set(gca, 'XScale', 'log')
        hold on

        % Add a line at the edge of the color patch
        plot3(myFreqs,repmat(1,size(myFreqs)),lgnYFitMatrix(whichStim,:),...
            '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.5);

        % Now the data
        lgnDataIndices = 1+(whichStim-1)*(nFreqs): ...
            (whichStim-1)*(nFreqs)+nFreqs;

        % Add the error bars
        yVals = lgnY(lgnDataIndices);
        ySEM = lgnYSEM(lgnDataIndices);
        for bb=1:length(yVals)
            plot3([studiedFreqs(bb) studiedFreqs(bb)],[1 1], ...
                [yVals(bb)-ySEM(bb), yVals(bb)+ySEM(bb)],...
                '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.0);
        end

        % Add the data symbols
        plot3(studiedFreqs,repmat(1,size(studiedFreqs)),yVals,...
            'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
            'MarkerSize',5,'MarkerEdgeColor','w','LineWidth',0.5);

        % Add a text label for the eccentricitiy
        text(150,1,0,'LGN',"HorizontalAlignment","left");

        % V1
        % Create a patch for the response
        X = [myFreqs fliplr(myFreqs)];
        Y = repmat(2,size(X));
        Z = [v1AvgFit, zeros(size(v1AvgFit))];
        p = fill3(X,Y,Z,plotColor{stimOrder(whichStim)});
        set(p,'edgecolor','none','facealpha',1);
        hold on

        % Add a line at the edge of the color patch
        plot3(myFreqs,repmat(2,size(myFreqs)),v1AvgFit,...
            '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.5);

        % Add the error bars
        for bb=1:length(v1YAvglow)
            plot3([studiedFreqs(bb) studiedFreqs(bb)],[2 2], ...
                [v1YAvglow(bb), v1YAvghigh(bb)],...
                '-','Color',lineColor{stimOrder(whichStim)},'LineWidth',1.0);
        end

        % Add the data symbols
        plot3(studiedFreqs,repmat(2,size(studiedFreqs)),v1AvgY,...
            'o','MarkerFaceColor',lineColor{stimOrder(whichStim)},...
            'MarkerSize',5,'MarkerEdgeColor','w','LineWidth',0.5);

        % Add a text label for the eccentricitiy
        text(150,2,0,'V1',"HorizontalAlignment","left");


    end

    % Clean up
    for ff=1:2
        figure(figHandles{ff});
        for ss=1:3
            nexttile(ss);
            view(0,viewAngles(ff));
            xlim([0.5 150])
            zlim([-1 7])
            a=gca;
            a.ZTick = [0,5];
            a.ZTickLabel = {'0','5'};
            a.XTick = [1,100];
            a.XTickLabel = {'1','100'};
            a.XTickLabelRotation = 0;
            a.XMinorTick = 'off';
            a.YAxis.Visible = 'off';
            box off
            if ss>1
                a.XAxis.Visible = 'off';
                a.ZAxis.Visible = 'off';
            end
        end
    end

    % Save the plots
    plotNamesPDF = {[subjects{whichSub} '_StackedModelFits.pdf' ],...
        [subjects{whichSub} '_AvgV1LGNStackedModelFits.pdf' ]};
    plotNamesPNG = {[subjects{whichSub} '_StackedModelFits.png' ],...
        [subjects{whichSub} '_AvgV1LGNStackedModelFits.png' ]};

    for ff=1:2
        figure(figHandles{ff});

        % Save just the vector elements
        for ss=1:3
            nexttile(ss);
            child_handles = allchild(gca);
            for cc=1:length(child_handles)
                if isgraphics(child_handles(cc),'patch')
                    set(child_handles(cc), 'Visible','off')
                else
                    set(child_handles(cc), 'Visible','on')
                end
            end
        end
        saveas(figHandles{ff},fullfile(savePath,plotNamesPDF{ff}));

        for ss=1:3
            nexttile(ss);

            child_handles = allchild(gca);
            for cc=1:length(child_handles)
                if ~isgraphics(child_handles(cc),'patch')
                    set(child_handles(cc), 'Visible','off')
                else
                    set(child_handles(cc), 'Visible','on')
                end
                a=gca;
                a.XAxis.Visible = 'off';
                a.YAxis.Visible = 'off';
                a.ZAxis.Visible = 'off';
            end
        end
        print(figHandles{ff},fullfile(savePath,plotNamesPNG{ff}),'-dpng','-r600');
    end
end
