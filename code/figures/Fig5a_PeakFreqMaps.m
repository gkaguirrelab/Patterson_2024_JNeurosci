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

% This is the threshold for the goodness of fit of the Watson model to the
% TTF in each vertex. We only display those voxels with this fVal or lower
fValThresh = 2;

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    if ss==1
        % Variables to hold the data across subjects
        vecByStimBySub = nan(nStims,nSubs,nVert);
        comboVec = nan(nSubs,nVert);
    end

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_WatsonFit_results.mat']);
    load(filePath,'fitResults')


    % Loop over stimulus directions and create a map of the peak frequency
    for whichStim = [3 1 2]

        % Find those vertices that had a positive response to this stimulus
        % direction
        fValSet = nan(size(results.R2));
        fValIdx = find(cellfun(@(x) ~isempty(x), fitResults.fVal));
        fValSet(fValIdx) = cellfun(@(x) x(whichStim), fitResults.fVal(fValIdx));

        % Identify those voxels with a positive response and an overall R2
        % of greater than the threshold, and an fVal below the threshold
        goodIdx = find(logical( (results.R2 > r2Thresh) .* (fValSet < fValThresh)  ));
        nGood = length(goodIdx);

        % Extract the peak amplitude and frequency for these vertices
        peakAmp = nan(nVert,1);
        peakAmp(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakAmp(goodIdx));
        peakFreq = nan(nVert,1);
        peakFreq(goodIdx) = cellfun(@(x) x(whichStim),fitResults.peakFreq(goodIdx));

        % Update the goodIdx
        goodIdx = find(~isnan(peakFreq));

        % save a peakAmp map
        newMap = templateImage;
        newMap.cdata = single(nan(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = single(peakAmp(goodIdx));
        newMap = ciftiMakePseudoHemi(newMap);
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakAmp.dtseries.nii']);
        cifti_write(newMap, fileOut);

        % save a peakFreq map
        newMap = templateImage;
        newMap.cdata = single(nan(size(fitResults.fVal)));
        newMap.cdata(goodIdx) = log10(single(peakFreq(goodIdx)));
        newMap = ciftiMakePseudoHemi(newMap);
        newMap.cdata(goodIdx) = 10.^(newMap.cdata(goodIdx));
        fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
        cifti_write(newMap, fileOut);

        vecByStimBySub(whichStim,ss,goodIdx) = log10(peakFreq(goodIdx));

        % Store a vector of Z-transformed peak-frequency values
        vec = log10(newMap.cdata(goodIdx));
        vec = vec - mean(vec);
        vec = vec ./ std(vec);
        comboVec(whichStim,goodIdx) = vec;

    end

    % save a combo peakFreq map
    newMap = templateImage;
    newMap.cdata = single(nan(size(fitResults.fVal)));
    newMap.cdata = single(mean(comboVec,1,'omitmissing'))';
    fileOut = fullfile(savePath,[subjectNames{ss} '_comboZpeakFreq.dtseries.nii']);
    cifti_write(newMap, fileOut);

    comboVecBySub(ss,:) = mean(comboVec,1,'omitmissing')';

end

% Write out an across-subject average peak-frequency map for each stimulus
for whichStim = [3 1 2]
    newMap = templateImage;
    newMap.cdata = single(nan(size(fitResults.fVal)));
    newVec = single(10.^mean(squeeze(vecByStimBySub(whichStim,:,:)),1,'omitmissing'))';
    newVec(isnan(newVec))=0;
    newMap.cdata = single(newVec);
    fileOut = fullfile(savePath,['AvgSubject_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
    cifti_write(newMap, fileOut);
end

% Create an across-subject average map
newMap = templateImage;
newMap.cdata = single(nan(size(fitResults.fVal)));
newMap.cdata = single(median(comboVecBySub,1,'omitmissing'))';
fileOut = fullfile(savePath,['AvgSubject_comboZpeakFreq.dtseries.nii']);
cifti_write(newMap, fileOut);

