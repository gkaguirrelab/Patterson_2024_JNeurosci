
%% Download the results file from the eventGain analysis

outputFileName = 'HEROgka1_eventGain_results.mat';

scratchSaveDir = getpref('forwardModelWrapper','flywheelScratchDir');
projectName = 'mtSinaiFlicker';
subjectName = 'vol2surfResults';
sessionName = 'vol2surfResults';


% Create the functional tmp save dir if it does not exist
saveDir = fullfile(scratchSaveDir,'v0','output');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Create a flywheel object
fw = flywheel.Flywheel(getpref('forwardModelWrapper','flywheelAPIKey'));

% The forwardModel results for one subject
searchStruct = struct(...
    'returnType', 'analysis', ...
    'filters', {{...
    struct('match', struct('project0x2elabel', projectName)), ...
    struct('match', struct('subject0x2ecode', subjectName)), ...
    struct('match', struct('session0x2elabel', sessionName)), ...
    struct('match', struct('analysis0x2elabel', 'forwardmodel')), ...
    }} ...
    );
analyses = fw.search(searchStruct);

% We should only find one analysis result for this search
if length(analyses)~=1
    error('Search failed to find a unique analysis')
end

% Get the analysis object
thisAnalysis = fw.getAnalysis(analyses{1}.analysis.id);

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
load(saveName,'results')


%% Fit the difference-of-exponentials model

nFreqs = 6;
freqs = [2 4 8 16 32 64];
freqsFit = logspace(log10(0.1),log10(100),1000);
nV = size(results.params,1);

% Anonymous function for a difference-of-exponentials model
doe = @(a,tau1,tau2,f) a.*(exp(-(f./tau1))-exp(-(f./tau2)));

% Variables to hold the results
fitPeakAmp = nan(nV,1);
fitPeakFreq = nan(nV,1);
fitR2 = nan(nV,1);

% Define the options for the search
options = optimset('fminsearch');
options.Display = 'off';

% Loop through the vertices / voxels
for vv = 1:nV
    if results.R2(vv)>0.1
        
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
        fitPeakAmp(vv) = a;
        fitPeakFreq(vv) = freqsFit(idx);        
        fitR2(vv) = corr(doe(x(1),x(2),x(3),freqs)',yVals');

    end
end


foo=1;


