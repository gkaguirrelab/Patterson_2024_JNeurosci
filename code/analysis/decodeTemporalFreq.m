function mriData = decodeTemporalFreq()

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

for ii=1:length(eccenDivs)-1
    eccenBins{ii}=[eccenDivs(ii),eccenDivs(ii+1)];
end

% The brain areas
areaLabels = {'v1','v2v3','lgn',};

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over areaLabels
    for aa = 1:length(areaLabels)

        switch areaLabels{aa}
            case 'lgn'
                areaIdx = LGNROI;
            case 'v1'
                areaIdx = (vArea==1);
            case 'v2v3'
                areaIdx = or(vArea==2,vArea==3);
        end

        inRangeIdx = areaIdx;

        % Filter the indices for goodness
        goodIdx = logical( (results.R2 > r2Thresh) .* inRangeIdx );

        % Loop over stimulus directions
        for dd = 1:length(stimulusDirections)

            % Set up a matrix to hold the voxel x acq results
            voxelByAcq = [];

            % Loop through the frequencies and obtain the set of values
            for ff = 1:nFreqs
                subString = sprintf(['f%dHz_' stimulusDirections{dd}],allFreqs(ff));
                idx = find(contains(stimLabels,subString));
                % Transpose to put voxels last
                voxelByAcq(ff,:,:) = results.params(goodIdx,idx)';
            end

            % Find all pairwise correlations between frequencies, getting
            % the average and SD across possible combinations
            pairSets = combinator(6,2,'c');


        end



    end

end


