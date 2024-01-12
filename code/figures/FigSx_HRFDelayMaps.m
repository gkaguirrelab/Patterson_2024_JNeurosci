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

% Get the flobs basis set with 100 ms resolution
deltaT = 0.1;
flobsbasis = returnFlobsVectors(deltaT);

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    if ss==1
        % Variables to hold the data across subjects
        vecBySub = nan(nSubs,nVert);
    end

    % Get the goodIdx
    goodIdx = find(results.R2 > r2Thresh);

    % Loop over the vertices and get the time to peak of the HRF
    for gg = 1:length(goodIdx)
        pp = results.params(goodIdx(gg),end-2:end);
        hrf = makeFlobsHRF(pp, flobsbasis);
        [~,idx] = max(hrf);
        timeToPeak{ss}(gg) = idx*deltaT;
    end

    vecBySub(ss,goodIdx) = timeToPeak{ss};

    % save the HRF delay map
    newMap = templateImage;
    newMap.cdata = single(nan(size(newMap.cdata)));
    newMap.cdata(goodIdx) = single(timeToPeak{ss});
    newMap = ciftiMakePseudoHemi(newMap);
    fileOut = fullfile(savePath,[subjectNames{ss} '_hrfDelayMap.dtseries.nii']);
    cifti_write(newMap, fileOut);

end


% save an across-subject HRF delay map
newMap = templateImage;
newMap.cdata = single(nan(size(results.fVal)));
newMap.cdata = single(mean(vecBySub,'omitnan'))';
newMap = ciftiMakePseudoHemi(newMap);
fileOut = fullfile(savePath,['AvgSubject_hrfDelayMap.dtseries.nii']);
cifti_write(newMap, fileOut);



%% Local functions



function hrf = makeFlobsHRF(x, flobsbasis)

% Create the HRF
hrf = flobsbasis*x';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

end