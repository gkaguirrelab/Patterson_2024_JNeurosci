function mriData = loadMRIResponseData()

%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
allFreqs = [0,2,4,8,16,32,64];
nFreqs = length(allFreqs);

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the LGN ROI.
tmpPath = fullfile(localDataDir,'retinoFiles','LGN_bilateral.dtseries.nii');
LGNROI = cifti_read(tmpPath); LGNROI = LGNROI.cdata;

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% The eccentricity bins that we will use to divide V1
eccenDivs = [0 90./(2.^(5:-1:0))];

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over stimulus directions
    for dd = 1:length(stimulusDirections)

        % Loop over eccentricity bins in V1 (plus LGN)
        for ee = 0:length(eccenDivs)

            % When ee==0, load the LGN data
            switch ee
                case 0 % The LGN
                    areaIdx = LGNROI;
                    eccFieldName = 'lgn';
                case length(eccenDivs) % Avg V1
                    % The eccentricity range for this bin
                    eccenRange = [eccenDivs(1) eccenDivs(end)];

                    % Find the vertices that we wish to analyze within V1
                    areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));

                    % The name we will use for this field in the data structure
                    eccFieldName = 'v1_avg';
                otherwise
                    % The eccentricity range for this bin
                    eccenRange = [eccenDivs(ee) eccenDivs(ee+1)];

                    % Find the vertices that we wish to analyze within V1
                    areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));

                    % The name we will use for this field in the data structure
                    eccFieldName = ['v1_ecc' num2str(ee)];
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


