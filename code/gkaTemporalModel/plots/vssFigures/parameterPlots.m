% scriptCreatePlots

% Housekeeping
clear

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% The identities of the stims and subjects
subjects = {'gka','asb'};
nSubs = length(subjects);
subLine = {'-',':'};
subSymbol = {'-',':'};
cellNames = {'midget','bistratified','parasol'};
lineColor={'r','b','k'};
nCells = length(cellNames);
paramNames = {'corner Freq','exponent','gain','rel. gain'};
nParams = 3;

% The eccentricities
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The colors used for the plots
faceAlpha = 0.4; % Transparency of the shaded error region

% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(800,400,'pt');
tiledlayout(1,4,'TileSpacing','tight','Padding','tight')

yLabels = {'Freq [Hz]','exponent','log gain','normed log gain'};
refCell = 3;

% Loop over params
for pp = 1:nParams+1
    nexttile(pp);

    for whichSub = 1:nSubs
        pMRI = storedSearchSeeds(subjects{whichSub});
        k = reshape(pMRI(2:end),nParams,nCells,nEccs);
        for cc = 1:nCells
            if pp>nParams
                vec = squeeze(k(pp-1,cc,:)) ./ squeeze(k(pp-1,refCell,:));
            else
                vec = squeeze(k(pp,cc,:));
            end
            plot(1:nEccs,vec,...
                ['o' subLine{whichSub}],'MarkerFaceColor',lineColor{cc},...
                'Color',lineColor{cc}, ...
                'MarkerSize',6,'MarkerEdgeColor','w','LineWidth',1);

            hold on
        end
    end
    title(paramNames{pp});
    xlabel('Eccentricity [deg]');
    ylabel(yLabels{pp})
    xlim([0 nEccs+1]);
    a=gca();
    a.XTick = 1:nEccs;
    a.XTickLabel = arrayfun(@num2str, round(studiedEccentricites), 'UniformOutput', 0);
    if pp >= 3
        a.YScale = 'log';
    end
end

% Save the plot
plotNamesPDF = 'paramPlots.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));


