clear
close all

% Place to save figures and load the results of the adapt analysis
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
nSubs = length(subjects);

figure

% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    fileName = fullfile(savePath,[subjectNames{ss} '_AvgV1_adaptMtSinai.mat']);
    load(fileName,'adaptResults')

        % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
    load(filePath,'results')

    % Get the adaptation params
    pp(ss,:) = adaptResults.params(end-11:end-3);

    plot(pp(ss,:),'*-')
    hold on

    % Report the R2 of the basic and adapt model
    fprintf([subjects{ss} ' - R2 basic, R2 adapt: %2.2f, %2.2f\n'],max(results.R2),max(adaptResults.R2));

end
