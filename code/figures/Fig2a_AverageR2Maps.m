clear
close all

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only display those voxels with this quality fit or better
r2Thresh = 0.1;

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
nSubs = length(subjects);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% Save a template map variable so we can create new maps below
tmpPath = fullfile(localDataDir,'MT.dtseries.nii');
templateImage = cifti_read(tmpPath);

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Get the vector of R2 values and nan values below the threshold
    vec = results.R2;
    vec(vec<r2Thresh) = nan;

    % save a map of the R2 fit value
    newMap = templateImage;
    newMap.cdata = single(vec);
    newMap = ciftiMakePseudoHemi(newMap);
    fileOut = fullfile(savePath,[subjectNames{ss} '_pseudoHemR2.dtseries.nii']);
    cifti_write(newMap, fileOut);

end

