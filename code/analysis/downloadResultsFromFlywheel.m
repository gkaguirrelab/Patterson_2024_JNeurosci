%% downloadResultsFromFlywheel

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};
analysisIDs = { {'60e9ea50a85492ed8f96cabd','60e9ea334ef89230db2b7021','60e9ea6dbceb4c0bc9e0767e'} , ...
    {'60e9eaa1a74445f40c56b123', '60e9ea85bd00f64426dd9301','60e9eabf4ef89230db2b7027'} };

% Create a location to save the retino files
retinoSaveDir = fullfile(localSaveDir,'retinoFiles');
mkdir(retinoSaveDir);

% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download and unzip the retino maps
retinoMapID = '5dc88aaee74aa3005e169380';
retinoFileName = 'TOME_3021_cifti_maps.zip';
tmpPath = fullfile(retinoSaveDir,retinoFileName);
fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
command = ['unzip -q -n '  escapeFileCharacters(tmpPath) ' -d ' escapeFileCharacters(retinoSaveDir)];
system(command);

% Create a location to save the results files
resultsSaveDir = fullfile(localSaveDir,'resultsFiles');
mkdir(resultsSaveDir);

% Loop through the subjects and directions and download results.mat
for ss = 1:2
        
    for dd = 1:3
        
        % Download the results file for this subject / direction
        fileStem = [subjectNames{ss} '_agtcOL_'];
        fileName = [fileStem 'results.mat'];
        tmpPath = fullfile(resultsSaveDir,[analysisLabels{dd} '_' fileName]);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{dd},fileName,tmpPath);
        
    end
end
        