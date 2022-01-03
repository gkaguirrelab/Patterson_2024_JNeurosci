%% downloadResultsFromFlywheel

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisIDs = { '6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75' };

% List of the output files to download
desiredOutputs = {'fig1.pdf','fig2.pdf','maps_cifti.zip','mtSinai_results.mat'};

% Create a location to save the retino files
retinoSaveDir = fullfile(localSaveDir,'retinoFiles');
mkdir(retinoSaveDir);

% Create a flywheel object. You need to set your flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download and unzip the retino maps
% retinoMapID = '5dc88aaee74aa3005e169380';
% retinoFileName = 'TOME_3021_cifti_maps.zip';
% tmpPath = fullfile(retinoSaveDir,retinoFileName);
% fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
% command = ['unzip -q -n '  escapeFileCharacters(tmpPath) ' -d ' escapeFileCharacters(retinoSaveDir)];
% system(command);

% Create a location to save the results files
resultsSaveDir = fullfile(localSaveDir,'resultsFiles');
mkdir(resultsSaveDir);

% Loop through the subjects and directions and download results.mat
for ss = 1:2

    fileStem = [subjectNames{ss} '_'];

    % Create a directory for this subject
    resultsSaveDir = fullfile(localSaveDir,[fileStem 'resultsFiles']);
    mkdir(resultsSaveDir);

    % Download the results file for this subject
    for ff = 1:length(desiredOutputs)
    fileName = [fileStem desiredOutputs{ff}];
    tmpPath = fullfile(resultsSaveDir,fileName);
    fw.downloadOutputFromAnalysis(analysisIDs{ss},fileName,tmpPath);
    end

end
