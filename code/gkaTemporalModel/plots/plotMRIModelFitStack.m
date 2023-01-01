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
spacing = 200;
stimOrder = [2 3 1];

plotColor={[0.85 0.85 0.85],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = [0,1,1];


for whichSub = 1:length(subjects)

patchHandles = {};

pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,myFreqs,rgcTemporalModel,paramCounts,modelType);

    figure('Renderer','painters');
    figuresize(600,300,'pt');

    for whichStim = 1:length(stimulusDirections)


        % Loop over eccentricities
        for ee=1:nEccs

            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % First plot the fit
            yVals = squeeze(v1YFitMatrix(whichStim,ee,:))';

            % Get the subplot and define the offset
            subplot(1,3,stimOrder(whichStim));
            set(gca, 'XScale', 'log')

            % Create a patch for the response
            X = [myFreqs fliplr(myFreqs)];
            Y = repmat(ee*spacing,size(X));
            Z = [yVals, zeros(size(yVals))];
            patchHandles{end+1} = fill3(X,Y,Z,plotColor{stimOrder(whichStim)});
            set(patchHandles{end},'edgecolor','none','facealpha',1);
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
                'MarkerSize',3,'MarkerEdgeColor','w','LineWidth',0.5);


            % Add a text label for the eccentricitiy
            text(150,ee*spacing,0,sprintf('%2.0f°',studiedEccentricites(ee)),"HorizontalAlignment","left");

        end


    end

    % Clean up
    for ss=1:3
        subplot(1,3,ss)
        view(0,70);
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
    end

    % Save the plot

    % Save just the vector elements
    for pp=1:length(patchHandles)
        set(patchHandles{pp}, 'Visible','off','HandleVisibility','off')
    end
    plotName = [subjects{whichSub} '_StackedModelFits.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

    for ss=1:3
        subplot(1,3,ss)
    child_handles = allchild(gca);
    for cc=1:length(child_handles)
        if ~isgraphics(child_handles(cc),'patch')
        set(child_handles(cc), 'Visible','off')
        else
        set(child_handles(cc), 'Visible','on')
        end
    end
    end
    plotName = [subjects{whichSub} '_StackedModelFits.png' ];
    print(fullfile(savePath,plotName),'-dpng','-r600');
    
end
