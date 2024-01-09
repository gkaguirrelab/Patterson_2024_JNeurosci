clear
close all

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
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
r2Thresh = 0.1;

figure

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
    lmsIdx = find(contains(stimLabels,'LMS'));

    % Get the z-score of the attention effect for these voxels
    pAtten=results.params(:,attenIdx);
    pLMS=results.params(:,lmsIdx);
    ztemp=mean(pAtten,2)./std(pAtten,[],2);
%    ztemp=mean(pLMS,2)./std(pLMS,[],2);
    ztemp(isnan(ztemp))=0;
    z(ss,:)=ztemp;

    % save the z attention map
    newMap = templateImage;
    newMap.cdata = single(squeeze(z(ss,:)))';
    newMap = ciftiMakePseudoHemi(newMap);
    fileOut = fullfile(savePath,[subjectNames{ss} '_attentionZmap.dtseries.nii']);
    cifti_write(newMap, fileOut);

    % Make a plot of the attention event effect as a function of
    % eccentricity. Load the V1 region and eccentricity for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    roiVol = cifti_read(filePath); roiVol = roiVol.cdata;
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_eccen.dscalar.nii']);
    eccenVol = cifti_read(filePath); eccenVol = eccenVol.cdata;

    % Get the good idx
    goodIdx = find(logical( (results.R2 > r2Thresh) .* (roiVol == 1)) .* (eccenVol > 0));

    % Scatter plot of eccentricity vs. z score of attention effect
    subplot(3,1,ss)
    x = log10(eccenVol(goodIdx));
    y = z(ss,goodIdx)';
    scatter(x,y,1,'k','.');
    hold on
    p = polyfit(x,y,1); 
    f = polyval(p,x); 
    plot(x,f,'-r') 
    xlim([0 2])
    xlabel('log10 Eccentricity');
    ylabel({'Attention response','[z score]'});
    title(sprintf([subjects{ss} ', slope = %2.2f'],p(1)));

end

% save an across-subject average attention map
newMap = templateImage;
newMap.cdata = single(zeros(size(results.fVal)));
z(z==0) = nan;
zAvg = mean(z,'omitnan');
zAvg(isnan(zAvg))=0;
newMap.cdata = single(zAvg)';
newMap = ciftiMakePseudoHemi(newMap);
fileOut = fullfile(savePath,['AvgSubject_attentionZmap.dtseries.nii']);
cifti_write(newMap, fileOut);
