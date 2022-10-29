function mriData = loadMRIResponseData()

%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
stimulusDirections = {'LminusM','S','LMS'};
allFreqs = [0,2,4,8,16,32,64];
nFreqs = length(allFreqs);

% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the subcortical ROIs
projectID = '5ca7803af546b60029ef118e';
subCorticalROIsFullNames = {'LGN_bilateral.dtseries.nii','thalamus_bilateral.dtseries.nii','midbrain_bilateral.dtseries.nii'};
subCorticalROIsLabels = {'LGN','thalamus','midbrain'};
for rr = 1:length(subCorticalROIsFullNames)
    tmpPath = fullfile(localSaveDir,'retinoFiles',subCorticalROIsFullNames{rr});
    fw.downloadFileFromProject(projectID,subCorticalROIsFullNames{rr},tmpPath);
    tmpRegion = cifti_read(tmpPath); tmpRegion = tmpRegion.cdata;
    str = [subCorticalROIsLabels{rr} 'ROI = tmpRegion;'];
    eval(str);
end

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% The eccentricity bins that we will use to divide V1
eccenDivs = [0 90./(2.^(5:-1:0))];

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localSaveDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over stimulus directions
    for dd = 1:length(stimulusDirections)

        % Loop over eccentricity bins in V1 (plus LGN)
        for ee = 0:length(eccenDivs)-1

            % When ee==0, load the LGN data
            if ee==0
                areaIdx = LGNROI;

                % The name we will use for this field in the data structure
                eccFieldName = 'lgn';
            else
                % The eccentricity range for this bin
                eccenRange = [eccenDivs(ee) eccenDivs(ee+1)];

                % Find the vertices that we wish to analyze within V1
                areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));

                % The name we will use for this field in the data structure
                eccFieldName = ['ecc' num2str(ee)];
            end

            % Get the indices that contain the desired results
            goodIdx = logical( (results.R2 > r2Thresh) .* areaIdx );

            % Loop through the frequencies and obtain the set of values
            rawVals = cell(1,nFreqs);
            for ff = 1:nFreqs
                subString = sprintf(['f%dHz_' stimulusDirections{dd}],allFreqs(ff));
                idx = find(contains(stimLabels,subString));
                rawVals{ff} = mean(results.params(goodIdx,idx),'omitnan');
            end

            % Adjust the values for the zero frequency
            adjustedVals = [];
            for ff = 2:nFreqs
                adjustedVals(:,ff-1) = rawVals{ff}-rawVals{1};
            end

            % Store these values in the data structure
            mriData.(shortNames{ss}).(stimulusDirections{dd}).(eccFieldName) = adjustedVals;

        end

    end
end


