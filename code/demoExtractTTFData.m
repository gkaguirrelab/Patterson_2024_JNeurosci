%% demoExtractTTFData
%
% This script demonstrates how to download the "results" files Flywheel and
% extract BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies


% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'LF','L-M','S'};
analysisIDs = { {'60e9ea6dbceb4c0bc9e0767e','60e9ea50a85492ed8f96cabd','60e9ea334ef89230db2b7021'} , ...
    {'60e9eabf4ef89230db2b7027', '60e9eaa1a74445f40c56b123', '60e9ea85bd00f64426dd9301'} };

% These are the frequencies that were studied
freqs = [2 4 8 16 32 64];
nFreqs = length(freqs);

% Set up a temporary location to save files we will need
scratchSaveDir = tempdir();

% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download and unzip the retino maps
retinoMapID = '5dc88aaee74aa3005e169380';
retinoFileName = 'TOME_3021_cifti_maps.zip';
saveDir = fullfile(scratchSaveDir,'v0','output');
mkdir(saveDir);
tmpPath = fullfile(saveDir,retinoFileName);
fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
command = ['unzip -q -n ' tmpPath ' -d ' saveDir];
system(command);

% Load the retino maps
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_varea.dtseries.nii'));
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_eccen.dtseries.nii'));
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_angle.dtseries.nii'));
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_sigma.dtseries.nii'));
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

% Download the results file for this subject / direction
fileStem = [subjectNames{ss} '_agtcOL_'];
fileName = [fileStem 'results.mat'];
tmpPath = fullfile(saveDir,[analysisLabels{dd} '_' fileName]);
fw.downloadOutputFromAnalysis(analysisIDs{ss}{dd},fileName,tmpPath);

% Load the result file into memory and delete the downloaded file
clear results
load(tmpPath,'results')

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


