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
plotColor = {'r','b','k'};
nCells = length(cellNames);
paramNames = {'corner Freq','exponent','gain'};
nParams = 3;

% Load the model fitting results
prefs.verbose = false;
projectDir = tbLocateProject('mriSinaiAnalysis',prefs);
resultFilePath = fullfile(projectDir,'data','cstResultsBootstrap.mat');
load(resultFilePath,'results')

% The eccentricities
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

% The colors used for the plots
faceAlpha = 0.1; % Transparency of the shaded error region

% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(800,400,'pt');
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

yLabels = {'Freq [Hz]','exponent','log gain'};
yLimSets = {[0 60],[0 2],10.^[1 4]};
refCell = 3;


for whichSub = 1:nSubs

    % Get the mean and the IQR of the parameters
    p = median(results.(subjects{whichSub}).p);
    pIQR = iqr(results.(subjects{whichSub}).p);
    pMat = sort(results.(subjects{whichSub}).p);
    pLow = p-pIQR/2;
    pHi = p+pIQR/2;

    yLimSets = {[0 70],[0 2.5],[0 200]};
    paramNames = {'corner Freq','exponent','gain'};
    subP = reshape(p(2:end),nParams,nCells,nEccs);
    subPlow = reshape(pLow(2:end),nParams,nCells,nEccs);
    subPhi = reshape(pHi(2:end),nParams,nCells,nEccs);

    for pp = 1:nParams

        % The gain values only have relative interpretations (not
        % absolute). We scale these here by 10 so that the log units are
        % all positive for the plot.
        if pp == 3
            subP = subP*10;
            subPlow = subPlow*10;
            subPhi = subPhi*10;
        end

        subplot(2,nParams,pp+3*(whichSub-1))

        for cc = 1:nCells

            thisVec = squeeze(subP(pp,cc,:))';
            thisLow = squeeze(subPlow(pp,cc,:))';
            thisHi = squeeze(subPhi(pp,cc,:))';

            % Add a patch for the error
            patch(...
                [log10(studiedEccentricites),fliplr(log10(studiedEccentricites))],...
                [ thisLow, fliplr(thisHi) ],...
                plotColor{cc},'EdgeColor','none','FaceColor',plotColor{cc},'FaceAlpha',faceAlpha);
            hold on
            plot(log10(studiedEccentricites),thisVec,[subLine{whichSub} plotColor{cc}],'LineWidth',2);
        end
        title(paramNames{pp});
        if pp == 2
            refline(0,1);
        end
        if pp == 3
            a = gca();
            a.YScale = 'log';
        end
        ylim(yLimSets{pp});
        ylabel(yLabels{pp});
    end

end

% Save the plot
plotNamesPDF = 'paramPlots.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));


