clear
close all

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

% This is the threshold for the goodness of fit of the Watson model to the
% TTF in each vertex. We only display those voxels with this fVal or lower
fValThresh = 2;


% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(400,400,'pt');
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

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
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValIdx = find(cellfun(@(x) ~isempty(x), fitResults.fVal));
        fValSet(fValIdx) = cellfun(@(x) x(whichStim), fitResults.fVal(fValIdx));
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));

        % Loop through the good vertices and calculate the total FI
        fisherInfo = [];
        for vv = 1:length(goodIdx)
            tuning = fitResults.yFitInterp{goodIdx(vv)};
            tuning = squeeze(tuning(whichStim,:));
            fi = ((diff(tuning).^2)./tuning(2:end));
            if isempty(fisherInfo)
                fisherInfo = fi;
            else
                fisherInfo = fisherInfo + fi;
            end
        end

        semilogx(interpFreqs(2:end),fisherInfo,'-','Color',stimPlotColors{whichStim});
        hold on

    end

end

%plotNamesPDF = 'ampAndFreqByVertexEccen.pdf';
%saveas(figHandle,fullfile(savePath,plotNamesPDF));

