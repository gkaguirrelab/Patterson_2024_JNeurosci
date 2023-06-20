% close all

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data');

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
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;


% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    figure

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

        % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
    load(filePath,'fitResults')

    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = 1:3

        % Find those vertices that had a positive response to this stimulus
        % direction
        posRespIdx = zeros(size(results.R2));
        posRespIdx(fitResults.eccDeg > 0) = cellfun(@(x) sum(x(whichStim,:))>1, fitResults.Y(fitResults.eccDeg > 0));

        % Identify those voxels with a positive response and an overall R2
        % of greater than the threshold
        goodIdx = find(logical( (results.R2 > r2Thresh) .* posRespIdx  ));
        nGood = length(goodIdx);

        % Extract the peak frequency for these vertices
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Filter the mis-fit vertices in which the max or minimum frequency
        % was identified
        peakFreq(peakFreq == 1) = nan;
        peakFreq(peakFreq > 40) = nan;

        % Update the goodIdx
        goodIdx = find(~isnan(peakFreq));

        % save a peakFreq map
        newMap = templateImage;
        newMap.cdata = single(zeros(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = single(peakFreq(goodIdx));
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
        cifti_write(newMap, fileOut);

        % Plot freq vs eccentricity
        x = log10(fitResults.eccDeg(goodIdx));
        x(x<0)=0;
        [x, sortedIdx] = sort(x);
        xq = 0:0.01:1.8;
        v = peakFreq(goodIdx);
        v = v(sortedIdx);
        plot(x,v,['.',stimPlotColors{whichStim}]);
        hold on
        sp = spaps(x,v,-100);
        vq = fnval(sp,xq);
        plot(xq,vq,['-' stimPlotColors{whichStim}],'LineWidth',3)


    end


end
