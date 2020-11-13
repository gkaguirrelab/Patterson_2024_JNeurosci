% Script that downloads the forwardModel results from the Mt Sinai flicker
% frequency experiments, and then fits a difference-of-exponentials model
% to the data. The resulting fits are saved back as maps.


% Save location for the maps
resultsSaveDir = '/Users/aguirre/Desktop/gka';
subjectName = 'HERO_gka';

% Analysis parameters
outputFileName = 'HEROgka1_eventGain_results.mat';
scratchSaveDir = getpref('forwardModelWrapper','flywheelScratchDir');
gearName = 'forwardmodel';
sessionIDs = {'5cab7166f546b6002af04970'};
analysisLabels = {'LightFlux','LMinusM','S'};

% Create the functional tmp save dir if it does not exist
saveDir = fullfile(scratchSaveDir,'v0','output');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Create a template image
%% REPLACE THIS WHEN WE HAVE THE FORWARD MODEL SAVE A TEMPLATE IMAGE
workbenchPath = getpref('forwardModelWrapper','wbCommand');
rawName = '/Users/aguirre/Downloads/task_sub-HEROgka1_ses-041416_task-LightFluxA_run-01_space-T1w_desc-preproc_bold_denoised_surfaces/ciftiFSLR_32k/task_sub-HEROgka1_ses-041416_task-LightFluxA_run-01_space-T1w_desc-preproc_bold_denoised_Atlas_s0.dtseries.nii';
templateImage = cifti_read(rawName, 'wbcmd', workbenchPath);
templateImage.cdata = templateImage.cdata(:,1);
   
% Create a flywheel object
fw = flywheel.Flywheel(getpref('forwardModelWrapper','flywheelAPIKey'));

% Get the analyses for this session
analyses = fw.getSessionAnalyses(sessionIDs{1});

% Loop through the analysis labels
for aa = 1:length(analysisLabels)
    
    % Find any analyses for this session that match the gearName and have a
    % label that starts with the analysisLabel string
    analysisIdx = and(...
        cellfun(@(x) startsWith(x.label,analysisLabels{aa}),analyses), ...
        cellfun(@(x) strcmp(x.gearInfo.name,gearName),analyses) );
    
    % Check that we have only one of these
    if sum(analysisIdx) > 1
        warning(['Found more than one forwardModel analyses with the label: ' analysisLabels{aa} '; skipping.' ]);
        continue
    elseif sum(analysisIdx) == 0
        warning(['Did not find a forwardModel analysis with the label: ' analysisLabels{aa} '; skipping.' ]);
        continue
    end
    
    analysisID = analyses{analysisIdx}.id;
    
    % Get the analysis object
    thisAnalysis = fw.getAnalysis(analysisID);
    
    % Find the file with the matching stem
    analysisFileMatchIdx = cellfun(@(x) strcmp(x.name,outputFileName),thisAnalysis.files);
    
    % Get some more information about the analysis and define a saveStem
    thisName = thisAnalysis.files{analysisFileMatchIdx}.name;
    saveName = fullfile(saveDir,thisName);
    
    % If the file has not already been downloaded, get it
    if ~exist(saveName,'file')
        % Inform the user
        fprintf(['Downloading: ' thisName '\n']);
        fprintf(['         to: ' saveDir '\n']);
        
        % Download the matching file to the rootSaveDir
        fw.downloadOutputFromAnalysis(thisAnalysis.id,thisName,saveName);
    end
    
    % Load the result file into memory
    clear results
    load(saveName,'results')
    
    % Delete the file on disk
    delete(saveName)

    % Fit the DoE model
    [results,fieldNames] = fitDoEModel(results);
        
    % Save the map results into images
    for ff = 1:length(fieldNames)
        
        % The initial, CIFTI space image
        outCIFTIFile = fullfile(resultsSaveDir, [subjectName '_' analysisLabels{aa} '_' fieldNames{ff} '.dtseries.nii']);
        outData = templateImage;
        outData.cdata = single(results.(fieldNames{ff}));
        outData.diminfo{1,2}.length = 1;
        cifti_write(outData, outCIFTIFile)
        
        % Create a pseudo-hemisphere version
        %ciftiMakePseudoHemi(outCIFTIFile, scratchSaveDir, 1, resultsSaveDir)

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
