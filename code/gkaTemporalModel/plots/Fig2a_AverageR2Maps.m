clear
close all

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
nSubs = length(subjects);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Save a template map variable so we can create new maps below
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
templateImage = cifti_read(tmpPath);

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % save a peakAmp map
    newMap = templateImage;
    newMap.cdata = results.R2;
    newMap = ciftiMakePseudoHemi(newMap);
    fileOut = fullfile(savePath,[subjectNames{ss} '_pseudoHemR2.dtseries.nii']);
    cifti_write(newMap, fileOut);

end

