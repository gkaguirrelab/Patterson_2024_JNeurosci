clear


% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
stimAlphas = [0.05 0.05 0.1];
nSubs = length(subjects);
nStims = length(stimulusDirections);

interpFreqs = logspace(log10(1),log10(100),501);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Save a template map variable so we can create new maps below
templateImage = cifti_read(tmpPath);

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.1;

plotOrder = [2,3,1];

% Loop through subjects and fit each vertex
for ss = 2:2%length(subjectNames)

    % Prepare the figures
    figHandle = figure('Renderer','painters');
    figuresize(600,200,'pt');
    tiledlayout(1,3,'Padding','tight')

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(savePath,[subjectNames{ss} '_WatsonFitUnconstrained_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = 1:3

        % Find those vertices that had a positive response to this stimulus
        % direction
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea > 0) .* (vArea < 4) ));

        maxVal = 0;

        % Loop through the good vertices and calculate the total FI
        nexttile(plotOrder(whichStim))
        fisherInfo = [];
        for vv = 1:length(goodIdx)
            tuning = fitResults.yFitInterp{goodIdx(vv)};
            tuning = squeeze(tuning(whichStim,:))*100;
            fi = ((diff(tuning).^2)./tuning(2:end));
            if isempty(fisherInfo)
                fisherInfo = fi;
            else
                fisherInfo = fisherInfo + fi;
            end
            if rand() < 0.05
                semilogx(interpFreqs(2:end),fi,'-','Color',[0.7,0.7,0.7],'LineWidth',0.4);
                maxVal = max([maxVal,max(fi)]);
                hold on
            end
        end
        ylim([0 0.3]);
        if plotOrder(whichStim) == 1
            ylabel('Fisher info')
        end
        a = gca();
        a.TickDir = 'out';
        box off
        a.XTickLabel = {'1','10','100'};

        if plotOrder(whichStim) > 1
            a.YTick = [];
        end

        yyaxis right
        semilogx(interpFreqs(2:end),fisherInfo,'-','Color',stimPlotColors{whichStim},'LineWidth',2);
        if plotOrder(whichStim) == 3
            ylabel('Total Fisher info')
        end
        if plotOrder(whichStim) < 3
            a.YTick = [];
        end

        ylim([0 300]);
        box off
        if plotOrder(whichStim) == 2; xlabel('Frequency [Hz]'); end

    end

end

plotNamesPDF = 'fisherInfoV1-3.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

