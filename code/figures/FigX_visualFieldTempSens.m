clear
close all

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
rangeValSet = {[5 25],[5 25],[5 25]};

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.25;

% This is the threshold for the goodness of fit of the Watson model to the
% TTF in each vertex. We only display those voxels with this fVal or lower
fValThresh = 1;

% Prepare the figures
figHandle = figure('Renderer','painters');
figuresize(400,600,'pt');
tiledlayout(nSubs,3,'TileSpacing','tight','Padding','tight')

% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

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

    % Initialize or load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFit_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValIdx = find(cellfun(@(x) ~isempty(x), fitResults.fVal));
        fValSet(fValIdx) = cellfun(@(x) x(whichStim), fitResults.fVal(fValIdx));
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh) .* (vArea == 1) ));

        % Assemble thge components for the field map plot. We assume that
        % the sigma size (in degrees of visual angle) of each vertex is
        % equivalent to the eccentricity multiplied by a scale factor of
        % 0.16, which is derived from Figure 7A of Kay et al 2013,
        % Compressive spatial summation in human visual cortex, J Neurophys
        vals = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));
        polarVals = polarMap(goodIdx);
        eccenVals = eccenMap(goodIdx);
        sigmaVals = eccenVals*0.16;
        rangeVals = rangeValSet{whichStim};

         nexttile;
        createFieldMap(vals,polarVals,eccenVals,sigmaVals,rangeVals);
        title([subjects{ss} ' - ' stimulusDirections{whichStim}]);

    end

end

plotNamesPDF = 'foo.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

