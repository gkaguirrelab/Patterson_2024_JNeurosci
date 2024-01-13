%% downloadResultsFromFlywheel

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};

% This forward model analysis used a smaller cortical mask and the LGN
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};

% This forwardModel analysis used the large cortical mask but no LGN
% analysisIDs = { '64a72e7ac71249914655c704','64a72ed7e54b4460ea4bf010' };


% List of the output files to download
desiredOutputs = {'fig1.pdf','fig2.pdf','maps_cifti.zip','mtSinai_results.mat'};

% Create a location to save the retino files
retinoSaveDir = fullfile(localDataDir,'retinoFiles');
mkdir(retinoSaveDir);

% Create a flywheel object. You need to set your flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download and unzip the retino maps. Currently commented out to avoid the
% time consuming re-download.
%{
    retinoMapID = '5dc88aaee74aa3005e169380';
    retinoFileName = 'TOME_3021_cifti_maps.zip';
    tmpPath = fullfile(retinoSaveDir,retinoFileName);
    fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
    command = ['unzip -q -n '  escapeFileCharacters(tmpPath) ' -d ' escapeFileCharacters(retinoSaveDir)];
    system(command);
%}

% Create a location to save the results files
resultsSaveDir = fullfile(localDataDir,'resultsFiles');
mkdir(resultsSaveDir);

% Loop through the subjects and directions and download results.mat
for ss = 1:2

    fileStem = [subjectNames{ss} '_'];

    % Create a directory for this subject
    resultsSaveDir = fullfile(localDataDir,[fileStem 'resultsFiles']);
    mkdir(resultsSaveDir);

    % Download the results file for this subject
    for ff = 1:length(desiredOutputs)
        fileName = [fileStem desiredOutputs{ff}];
        tmpPath = fullfile(resultsSaveDir,fileName);
        fw.downloadOutputFromAnalysis(analysisIDs{ss},fileName,tmpPath);
    end

end
