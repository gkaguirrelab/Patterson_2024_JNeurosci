

%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.05;
nBoots = 50;

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
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

            % Store the confusion Matrix
            avConfuse(dd,ee,ss,:,:) = confusionMatrix;

        end % directions

    end % Eccentricity Bins

end % Subjects


% Get the achromatic and chromatic confusion matrices
achrConfuse = squeeze(avConfuse(3,:,:,:,:));
chromConfuse = squeeze(mean(avConfuse(1:2,:,:,:,:),1));

% Plot the correct decoding performance for 64 Hz stimuli across
% eccentricity
% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
subMarkers = {'^','square','o'};
subMarkerSize = [9,11,8];
subLines = {'-','--',':'};
setDirections = {'achromatic','chromatic'};

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.85]};
lineColor={'k','m'};
dirPlotShift = 0.1;

% Prepare the figure
figHandle = figure('Renderer','painters');
figuresize(400,200,'pt');

for ss = 1:nSubs
    x = 1:length(eccenBins);
    x = x + (ss-2)*dirPlotShift;
    y = squeeze(achrConfuse(:,ss,6,6));
    plot(x,y,subMarkers{ss},...
        'MarkerEdgeColor','none','MarkerFaceColor',plotColor{1},...
        'MarkerSize',subMarkerSize(ss))
    hold on
    y = squeeze(chromConfuse(:,ss,6,6));
    plot(x,y,subMarkers{ss},...
        'MarkerEdgeColor','none','MarkerFaceColor',plotColor{2},...
        'MarkerSize',subMarkerSize(ss))
end

x = 1:length(eccenBins);
% Median across subjects
y = median(squeeze(achrConfuse(:,:,6,6)),2);
plot(x,y,'-','Color',lineColor{1},'LineWidth',2);
plot(x,y,'o',...
    'MarkerFaceColor','none','MarkerEdgeColor',lineColor{1},...
    'MarkerSize',15,'LineWidth',2)
y = median(squeeze(chromConfuse(:,:,6,6)),2);
plot(x,y,'-','Color',lineColor{2},'LineWidth',2);
plot(x,y,'o',...
    'MarkerFaceColor','none','MarkerEdgeColor',lineColor{2},...
    'MarkerSize',15,'LineWidth',2)

a = gca();
a.XTick = x;
a.XTickLabel = {'2','4','8','16','32','64'};
xlabel('Eccentricity bin center [degrees]');
ylabel({'Proportion correct','decoding 64 Hz flicker'})
ylim([0,1.1]);
a.YTick = [0,0.5,1];
plot([0.5 6.5],[1/6 1/6],':k');
box off
drawnow


% Save the plot
plotNamesPDF = 'Fig7-2_64HzDecodeByEccen.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

foo=1;
