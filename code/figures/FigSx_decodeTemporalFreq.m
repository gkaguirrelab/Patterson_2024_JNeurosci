

%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.05;
nBoots = 50;

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
subMarkers = {'^','square','o'};
subMarkerSize = [9,11,8];
subLines = {'-','--',':'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(allFreqs);

% The eccentricity bins that we will use to divide V1
eccenDivs = [0 90./(2.^(5:-1:0))];
for ii=1:length(eccenDivs)-1
    eccenBins{ii}=[eccenDivs(ii),eccenDivs(ii+1)];
end

% The brain areas
areaLabels = {'v1'};

% Color map
cmap = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';

% Set the column order in the plot
colOrder = [2 3 1];

% Loop through subjects
for ss = 1:nSubs

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Load the ROI files for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    vArea = cifti_read(filePath); vArea = vArea.cdata;
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_eccen.dscalar.nii']);
    eccenVol = cifti_read(filePath); eccenVol = eccenVol.cdata;

    % Loop over areaLabels
    for ee = 1:length(eccenBins)

        figHandle = figure('Name',[subjects{ss} ' - eccenBin: ' num2str(ee)]);
        t = tiledlayout(2,3);
        t.TileSpacing = 'compact';
        t.Padding = 'compact';

        % Filter the indices for goodness
        goodIdx = logical( (results.R2 > r2Thresh) .* (vArea==1) .* (eccenVol>=eccenBins{ee}(1)) .* (eccenVol<=eccenBins{ee}(2)) );

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
            confusionMatrix = zeros(6,6);
            trainingTargets = [];
            for pp = 1:size(pairSets,1)
                setA = pairSets(pp,:);
                setB = idxVec(~ismember(idxVec,setA));
                % Loop through the training set and build the training
                % matrix
                for ii=1:6
                    trainingTargets(ii,:) = mean(squeeze(voxelByAcq(dd,ii,setA,:)))';
                end
                for ii=1:6
                    for jj=1:6
                        vecA = mean(squeeze(voxelByAcq(dd,ii,setA,:)))';
                        vecB = mean(squeeze(voxelByAcq(dd,jj,setB,:)))';
                        corrMat(pp,ii,jj) = corr(vecA,vecB);
                        % See which training target is most similar to vecB
                        [~,matchIdx] = max(corr(trainingTargets',vecB));
                        confusionMatrix(jj,matchIdx) = confusionMatrix(jj,matchIdx)+1;
                    end
                end
            end
            confusionMatrix = confusionMatrix ./ sum(confusionMatrix,2);

            % Plot the correlation matrix
            nexttile(colOrder(dd));
            im = squeeze(mean(corrMat));
            im = round(im * 128 + 128);
            image(im);
            colormap(cmap)

            axis square
            if colOrder(dd) == 1
                plotCleanUp(allFreqs);
            else
                axis off
            end
            stdMap = squeeze(std(corrMat));
            meanSEM = mean(stdMap(:));
            title(['Pearson - ' stimulusDirections{dd} sprintf(' sem=%2.2f',meanSEM)]);
            drawnow

            % Plot the confusion matrix
            nexttile(colOrder(dd)+3);
            im = confusionMatrix;
            im = round(im * 128 + 128);
            image(im);
            colormap(cmap)

            axis square
            axis off
            title(['confusion - ' stimulusDirections{dd}]);
            drawnow

        end % directions

        nexttile(6)
        cb = colorbar('southoutside');
        cb.Ticks = linspace(1,256,5);
        cb.TickLabels=arrayfun(@(x) {num2str(x)},[-1 -0.5 0 0.5 1]);
        cb.TickLength = [0 0];
        box off

    end % Areas

end % Subjects



%% LOCAL FUNCTION
function plotCleanUp(allFreqs)
a = gca();
a.XTick = 1:length(allFreqs)-1;
a.YTick = 1:length(allFreqs)-1;
a.XTickLabels = arrayfun(@(x) {num2str(x)},allFreqs(2:end));
a.YTickLabels = arrayfun(@(x) {num2str(x)},allFreqs(2:end));
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
xlabel('freq [Hz]');
ylabel('freq [Hz]');
end
