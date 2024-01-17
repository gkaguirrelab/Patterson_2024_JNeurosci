
clear
close all

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names and stimulus directions
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
shortNames = {'gka','asb','cgp'};
stimClassSet = {...
    {'f2Hz_LminusM','f2Hz_S'},...
    {'f4Hz_LminusM','f4Hz_S'},...
    {'f8Hz_LminusM','f8Hz_S'},...
    {'f16Hz_LminusM','f16Hz_S'},...
    {'f32Hz_LminusM','f32Hz_S'},...
    {'f64Hz_LminusM','f64Hz_S'},...
    {'f2Hz_LMS'},...
    {'f4Hz_LMS'},...
    {'f8Hz_LMS'},...
    {'f16Hz_LMS'},...
    {'f32Hz_LMS'},...
    {'f64Hz_LMS'};
    };

% Loop through the subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
    load(filePath,'results')

    % Define model components
    datats = results.data.datats';
    for ii = 1:36
        data{ii}=datats((ii-1)*336+1:ii*336);
    end
    stimulus = results.model.inputs{2};
    modelOpts = results.model.opts;
    tr = results.meta.tr;
    modelClass = 'mtSinaiAdapt';

    % Special case of handling the stimtime variable for cgp
    if ss==3
        modelOpts = [modelOpts(1:9) {modelOpts(10:end)}];
    end

    % Add the stimClassSet we wish to model
    modelOpts{end+1}='stimClassSet';
    modelOpts{end+1}=stimClassSet;

    % Call the forwardModel
    adaptResults = forwardModel(data,stimulus,tr,...
        'averageVoxels',true,...
        'verbose',false,...
        'modelClass',modelClass,...
        'modelOpts',modelOpts);

    % Save the results
    fileName = fullfile(savePath,[subjectNames{ss}  '_AvgV1_adaptMtSinai.mat']);
    save(fileName,'adaptResults');

end

