function decodeTemporalFreq()

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
areaLabels = {'v1','v2v3'};

% Loop through subjects
for ss = 1:2

    figure('Name',shortNames{ss})
    t = tiledlayout(2,3);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over areaLabels
    for aa = 1:2

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
            temp = [];
            for ff = 1:nFreqs
                subString = sprintf(['f%dHz_' stimulusDirections{dd}],allFreqs(ff));
                idx = find(contains(stimLabels,subString));
                % Transpose to put voxels last
                paramSet = results.params(goodIdx,idx)';
                % Subtract the zero frequency condition
                temp(ff,:,:) = paramSet;
            end

            % Subtract off the zero hz condition
            temp = temp - temp(1,:,:);

            % Store the non-zero conditiomns
            voxelByAcq(dd,:,:,:) = temp(2:end,:,:);

            % Loop through the repetitions and standardize the across
            % frequency values for each voxel
            for rr = 1:size(voxelByAcq,3)
                myMat = squeeze(voxelByAcq(dd,:,rr,:));
                myMat = myMat - mean(myMat);
                myMat = myMat ./ std(myMat);
                % Mean-center the values across
                % voxels at each frequency
                myMat = myMat - mean(myMat,2);
                voxelByAcq(dd,:,rr,:) = myMat;
            end


            % Find all pairwise correlations between frequencies, getting
            % the average and SD across possible combinations of
            % acquisitions
            pairSets = combinator(size(voxelByAcq,3),size(voxelByAcq,3)/2,'c');
            idxVec = 1:size(voxelByAcq,3);
            corrMat = [];
            for pp = 1:size(pairSets,1)
                for ii=1:6
                    for jj=1:6
                        setA = pairSets(pp,:);
                        setB = idxVec(~ismember(idxVec,setA));
                        vecA = mean(squeeze(voxelByAcq(dd,ii,setA,:)))';
                        vecB = mean(squeeze(voxelByAcq(dd,jj,setB,:)))';
                        corrMat(pp,ii,jj) = corr(vecA,vecB);
                    end
                end

            end

            % Plot the correlation matrix
            nexttile;
            im = squeeze(mean(corrMat));
            im = round(im * 128 + 128);
            image(im);
            colormap(redblue);

            axis square
            axis off
            title([areaLabels{aa} '-' stimulusDirections{dd}]);
            drawnow

        end % directions

    end % Areas


end % Subjects


end % function


