clear
close all

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

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

% Save a template map variable so we can create new maps below
tmpPath = fullfile(localDataDir,'MT.dtseries.nii');
templateImage = cifti_read(tmpPath);

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.0;

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Get the goodIdx
    goodIdx = ~isnan(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Get the indices of the covariates that model attention events
    attenIdx = find(contains(stimLabels,'attention'));

    % Get the z-score of the attention effect for these voxels
    p=results.params(goodIdx,attenIdx);
    z=mean(p,2)./std(p,[],2);

    % save the z attention map
    newMap = templateImage;
    newMap.cdata = single(zeros(size(results.fVal)));
    newMap.cdata(goodIdx) = single(z);
    newMap = ciftiMakePseudoHemi(newMap);
    fileOut = fullfile(savePath,[subjectNames{ss} '_attentionZmap.dtseries.nii']);
    cifti_write(newMap, fileOut);
end
