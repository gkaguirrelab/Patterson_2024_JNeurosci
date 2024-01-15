
clear
close all

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names and stimulus directions
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
shortNames = {'gka','asb','cgp'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};
stimClassSetSource = {'f2Hz_XXX','f4Hz_XXX','f8Hz_XXX','f16Hz_XXX','f32Hz_XXX','f64Hz_XXX'};

% Loop through the subjects
for whichStim = 1:length(directions)
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
        stimClassSet = strrep(stimClassSetSource,'XXX',directions{whichStim});
        modelOpts{end+1}='stimClassSet';
        modelOpts{end+1}=stimClassSet;

        % Call the forwardModel
        adaptResults = forwardModel(data,stimulus,tr,...
            'averageVoxels',true,...
            'verbose',false,...
            'modelClass',modelClass,...
            'modelOpts',modelOpts);

        % Save the results
        fileName = fullfile(savePath,[subjectNames{ss} '_' directions{whichStim} '_AvgV1_adaptMtSinai.mat']);
        save(fileName,'adaptResults');

    end
end


%% Local Functions
function f=portraitFigure()
f = figure();
set(f,...
    'PaperPosition',[0 0 11 8.5],...
    'PaperSize',[11 8.5000],...
    'PaperOrientation','portrait',...
    'Position',[543 336 791 611],...
    'OuterPosition',[543 336 791 690],...
    'InnerPosition',[543 336 791 611]);

end

