% Script that downloads the forwardModel results from the Mt Sinai flicker
% frequency experiments, and then fits a difference-of-exponentials model
% to the data. The resulting fits are saved back as maps.


% Save location for the maps
subjectNames = {'HEROgka1','HEROasb1'};
analysisIDs = { {'602ae04bf32ab1cf2fe2d7a5','602ae042ae50da7a462c12d6','602ae03afae38976c793401b'} , ...
    {'60313a5e2471d261bd0b62e8', '60313a5593240e7ae82c122b', '60313a4cd184138a20933f3c'} };

% Analysis parameters
scratchSaveDir = getpref('forwardModelWrapper','flywheelScratchDir');
analysisLabels = {'LF','L-M','S'};

% Create the functional tmp save dir if it does not exist
saveDir = fullfile(scratchSaveDir,'v0','output');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Create a flywheel object
fw = flywheel.Flywheel(getpref('forwardModelWrapper','flywheelAPIKey'));

% Loop over subjects
for ss = 1: length(subjectNames)
    
    % Set up the paths for this subject
    fileStem = [subjectNames{ss} '_eventGain_'];    
    resultsSaveDir = ['/Users/aguirre/Desktop/' subjectNames{ss}];
    mkdir(resultsSaveDir);
    
    % Loop over analysis IDs
    for aa = 1:length(analysisIDs{ss})
        
        % Download the results file
        fileName = [fileStem 'results.mat'];
        tmpPath = fullfile(saveDir,fileName);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{aa},fileName,tmpPath);
        
        % Load the result file into memory and delete the downloaded file
        clear results
        load(tmpPath,'results')
        delete(tmpPath)
        
        % Download the templateImage file
        fileName = [fileStem 'templateImage.mat'];
        tmpPath = fullfile(saveDir,fileName);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{aa},fileName,tmpPath);
        
        % Load the result file into memory and delete the downloaded file
        clear templateImage
        load(tmpPath,'templateImage')
        delete(tmpPath)
        
        % Fit the DoE model
        [results,fieldNames] = fitDoEModel(results);
        
        % Save the map results into images
        for ff = 1:length(fieldNames)
            
            % The initial, CIFTI space image
            outCIFTIFile = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} '_' fieldNames{ff} '.dtseries.nii']);
            outData = templateImage;
            outData.cdata = single(results.(fieldNames{ff}));
            outData.diminfo{1,2}.length = 1;
            cifti_write(outData, outCIFTIFile)
            
        end
        
    end
    
end


function [results,fieldNames] = fitDoEModel(results)

%% Fit the difference-of-exponentials model
nFreqs = 6;
freqs = [2 4 8 16 32 64];
freqsFit = logspace(log10(0.1),log10(100),1000);
nV = size(results.params,1);

% Anonymous function for a difference-of-exponentials model
doe = @(a,tau1,tau2,f) a.*(exp(-(f./tau1))-exp(-(f./tau2)));

% Variables to hold the results
fitPeakAmpDoE = nan(nV,1);
fitPeakFreqDoE = nan(nV,1);
fitR2DoE = nan(nV,1);

% Define the options for the search
options = optimset('fminsearch');
options.Display = 'off';

% Loop through the vertices / voxels
for vv = 1:nV
    if results.R2(vv)>0.05
        
        % Get the beta values
        yVals = results.params(vv,1:nFreqs);
        
        % Handle a negative offset
        if min(yVals)<0
            offset = min(yVals);
            yVals = yVals-offset;
        else
            offset = 0;
        end
        
        % Set up the objective
        myObj = @(x) norm(yVals - doe(x(1),x(2),x(3),freqs));
        
        % Fit
        x = fminsearch(myObj,[20,1,2],options);
        myFit = doe(x(1),x(2),x(3),freqsFit)+offset;
        
        [a,idx] = max(myFit);
        fitPeakAmpDoE(vv) = a;
        fitPeakFreqDoE(vv) = freqsFit(idx);
        fitR2DoE(vv) = corr(doe(x(1),x(2),x(3),freqs)',(yVals+offset)');
        
    end
end

% Place the analysis in the results variable
results.fitPeakAmpDoE = fitPeakAmpDoE;
results.fitPeakFreqDoE = fitPeakFreqDoE;
results.fitR2DoE = fitR2DoE;
fieldNames = {'fitPeakAmpDoE','fitPeakFreqDoE','fitR2DoE'};

end
