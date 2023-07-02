clear
close all

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
stimAlphas = [0.05 0.05 0.1];
nSubs = length(subjects);
nStims = length(stimulusDirections);

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

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValSet(fitResults.eccDeg > 0) = cellfun(@(x) x(whichStim), fitResults.fVal(fitResults.eccDeg > 0));

        % Identify those voxels with a positive response and an overall R2
        % of greater than the threshold, and an fVal below the threshold
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh)  ));
        nGood = length(goodIdx);

        % Extract the peak amplitude and frequency for these vertices
        peakAmp = nan(nVert,1);
        peakAmp(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakAmp(goodIdx));
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Update the goodIdx
        goodIdx = find(~isnan(peakFreq));

        % save a peakAmp map
        newMap = templateImage;
        newMap.cdata = single(zeros(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = single(peakAmp(goodIdx));
        newMap = ciftiMakePseudoHemi(newMap);
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakAmp.dtseries.nii']);
        cifti_write(newMap, fileOut);

        % save a peakFreq map
        newMap = templateImage;
        newMap.cdata = single(zeros(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = single(peakFreq(goodIdx));
        newMap = ciftiMakePseudoHemi(newMap);
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
        cifti_write(newMap, fileOut);

    end

end

