% scriptCreatePlots

% Housekeeping
clear

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% The identities of the stims and subjects
subjects = {'gka','asb'};
nSubs = length(subjects);
subSymbol = {'o','s'};
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
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

yLabels = {'Freq [Hz]','exponent','log gain'};
yLimSets = {[0 60],[0 2],10.^[0 3]};
refCell = 3;

% Loop over params
for pp = 1:nParams
    nexttile(pp);

    for whichSub = 1:nSubs

        % Get the LGN parameters
        lgnResultFile = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'cstResultsLGN.mat');
        load(lgnResultFile,'results');
        thisP = results.(subjects{whichSub}).p;

        for cc = 1:nCells
            val = thisP(1+(cc-1)*3+pp);

            % Given that the gain parameter has arbitrary units, we scale
            % the values so that they range from 1 --> larger.
            if pp==3
                val = val * 15;
            end

            plot(1,val,...
                subSymbol{whichSub},'MarkerFaceColor',lineColor{cc},...
                'Color',lineColor{cc}, ...
                'MarkerSize',15,'MarkerEdgeColor','w');
            hold on

        end
       
    end
    title(paramNames{pp});
    ylabel(yLabels{pp})
    xlim([0 nEccs+1]);
    a=gca();
    if pp == 2
        plot([0.5 1.5],[1 1],'-k');
    end
    if pp >= 3
        a.YScale = 'log';
    end
    xlim([0.75 1.25])
    ylim(yLimSets{pp});
end

% Save the plot
plotNamesPDF = 'paramPlots_wLGN.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));


