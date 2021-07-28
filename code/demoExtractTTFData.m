%% demoExtractTTFData
%
% This script demonstrates how to download the "results" files Flywheel and
% extract BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'LF','L-M','S'};

% These are the frequencies that were studied
freqs = [2 4 8 16 32 64];
nFreqs = length(freqs);

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Create a "subcortical" map
subcorticalMap = zeros(size(vArea));
subcorticalMap(1:26298)=1;

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.25;

% This is the visual area and eccentricity range to grab. The visual areas
% are: V1 = 1, V2 = 2, V3 = 3, hV4/LO = [4 5], MT/MST = [8 9]
area = 1;
eccenRange = [0 20];

% Which subject?
ss = 2;

% Which stimulus direction?
dd = 1;

% Load the results file for this subject / direction
filePath = fullfile(localSaveDir,'resultsFiles',[analysisLabels{dd} '_' subjectNames{ss} '_agtcOL_results.mat']);
load(filePath,'results')

% Find the vertices that we wish to analyze
goodIdx = logical( (results.R2 > r2Thresh) .* (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2)) );

% Get the beta values for these indices
yVals = results.params(goodIdx,1:nFreqs+1);

% Get the mean betas across vertices / voxels
yValsMean = nanmean(yVals);

% Obtain the set of beta valyes, relative to the baseline condition
y = yValsMean(2:end)-yValsMean(1);

% Plot the data
figure
semilogx(freqs,y,'-k');
hold on
semilogx(freqs,y,'*r');

% Clean up the plot
title(sprintf([shortNames{ss} '-' analysisLabels{dd} ' v = %d'],sum(goodIdx)))
ylabel('BOLD response [% change]')
xlabel('Frequency [hz]');


