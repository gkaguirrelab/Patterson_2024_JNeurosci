clear

%% NEED TO FIX THIS BY LOOPING THROUGH VERTICES AND GENERATING THE YFIT INTERP FOR THE WATSON FIT PARAMS

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
stimAlphas = [0.05 0.05 0.1];
nSubs = length(subjects);
nStims = length(stimulusDirections);
studiedFreqs = [2 4 8 16 32 64];

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.1;

plotOrder = [2,3,1];

% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Prepare the figures
    figHandle = figure('Renderer','painters');
    figuresize(400,200,'pt');
    tiledlayout(1,3,'TileSpacing','compact','Padding','tight')

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Load the retino maps for this subject
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    vArea = cifti_read(tmpPath); vArea = vArea.cdata;
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_eccen.dscalar.nii']);
    eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_angle.dscalar.nii']);
    polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFit_results.mat']);
    load(filePath,'fitResults')

    lowFreqIdx = 77;
    hiFreqIdx = 451;
    interpFreqs = logspace(log10(1),log10(100),501);
    interpFreqs=interpFreqs(77:451);

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = 1:3

        % Find those vertices that had a positive response to this stimulus
        % direction
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea > 0) .* (vArea < 2) ));

        maxVal = 0;

        % Loop through the good vertices and calculate the total FI
        nexttile(plotOrder(whichStim))
        fisherInfo = [];
        for vv = 1:length(goodIdx)
            signal = fitResults.yFitInterp{goodIdx(vv)};
            noise = fitResults.stdFitInterp{goodIdx(vv)};
            signal = squeeze(signal(whichStim,:))*100;
            noise = squeeze(noise(whichStim,:))*100;
            tuning = signal ./ noise;
            tuning = tuning(77:451);
            fi = tuning(2:end);
            %            fi = ((diff(tuning).^2)./tuning(2:end));
            if isempty(fisherInfo)
                fisherInfo = fi;
            else
                fisherInfo = fisherInfo + fi;
            end
            if rand() < 0.05
                semilogx(interpFreqs(2:end),log10(fi),'-','Color',[0.7,0.7,0.7],'LineWidth',0.4);
                maxVal = max([maxVal,max(fi)]);
                hold on
            end
        end
        ylim([-6 2]);
        if plotOrder(whichStim) == 1
            ylabel('log Fisher information')
        end
        a = gca();
        a.TickDir = 'out';
        a.XTick = studiedFreqs;
        xlim([2 64]);
        box off

        if plotOrder(whichStim) > 1
            a.YTick = [];
        end

        yyaxis right
        semilogx(interpFreqs(2:end),log10(fisherInfo),'-','Color',stimPlotColors{whichStim},'LineWidth',3);
        if plotOrder(whichStim) == 3
            ylabel('log total Fisher information')
        end

        [~, maxIdx] = max(fisherInfo);
        vec = fisherInfo; vec(1:maxIdx) = nan;
        fiElbow = log10(max(fisherInfo)) - 0.2;
        [~,elbowIdx] = find(log10(vec) < fiElbow,1,"first");
        x = interpFreqs(elbowIdx+1);
        y = log10(fisherInfo(elbowIdx));
        semilogx(x,y,'|','Color',stimPlotColors{whichStim},'MarkerSize',15);
        text(10.^(log10(x)-0.1),3.75,sprintf('%d Hz',round(x)));

        a.YTick = 1:4;
        if plotOrder(whichStim) < 3
            a.YTick = [];
        end

        ylim([1 4]);
        box off
        if plotOrder(whichStim) == 2; xlabel('Frequency [Hz]'); end

    end

    plotNamesPDF = [subjectNames{ss} '_fisherInfoV1.pdf'];
    saveas(figHandle,fullfile(savePath,plotNamesPDF));

end


