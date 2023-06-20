% Validate the fit for vertex 53477 for subject GKA


%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;
chunkSize = 100;

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% We try and improve the temporal model fit for fVal results above this
refitThresh = 5.0;

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Fixed features of the model
nAcqs = 12; nCells = 3; nParams = 3;

% The "quality" parameter of the low-pass filter. In initial analyses we
% found that using a quality ("Q") parameter of unity for the low-pass
% filter stage produces a good result for both subjects, so we lock that in
% analysis here across subjects and eccentricities.
Q = 1.0;

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(studiedFreqs);
interpFreqs = logspace(log10(2),log10(64),11);
freqsForPlotting = logspace(0,2,50);

% BADs Options. Indicate that the objective function is deterministic, and
% handle verbosity
optionsBADS.UncertaintyHandling = 0;
if verbose
    optionsBADS.Display = 'iter';
else
    optionsBADS.Display = 'off';
end

% The optimization toolbox is currently not available for Matlab running
% under Apple silicon. Detect this case and tell BADS so that it doesn't
% issue a warning
V = ver;
if ~any(strcmp({V.Name}, 'Optimization Toolbox'))
    optionsBADS.OptimToolbox = 0;
end

% A function to vectorize
vectorize = @(x) x(:);

% Save the warning state
warnstate = warning();

% Set the p0
p0 = [Q 40 2 0.1 20 1.0 10 30 1.0 1];
p1 = [Q 20 1 1.0 10 0.5 10 10 0.5 1];

% Bounds on Q, corner frequency, exponentiation, gain.
lb = [Q 1 0.1 0.01 1 0.1 0.01 1 0.1 0.01];
ub = [Q 100 3.0 100 100 3.0 100 100 3.0 100];


%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the LGN ROI.
tmpPath = fullfile(localDataDir,'retinoFiles','LGN_bilateral.dtseries.nii');
LGNROI = cifti_read(tmpPath); LGNROI = LGNROI.cdata;


%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % How many vertices total?
    nVert = length(results.fVal);

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Fit any voxel with an R2 fit above the threshold
    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea==1) ));
    nGood = length(goodIdx);

    % Initialize or load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
    if isfile(filePath)
        load(filePath,'fitResults')
    else
        fitResults = [];
        fitResults.p = nan(nVert,length(p0));
        fitResults.fVal = nan(nVert,1);
        fitResults.Y = cell(nVert,1);
        fitResults.yFit = cell(nVert,1);
        fitResults.eccDeg = nan(nVert,1);
        fitResults.polarAngle = nan(nVert,1);
        fitResults.area = nan(nVert,1);
    end

    % We chunk the parpool search. This is because we need to occasionally
    % shut down the par pool and re-open it due to memory leaks from the
    % symbolic toolbox
    nChunks = ceil(nGood/chunkSize);

    for cc = 3:nChunks

        % Set the start and end points of the goodIdx we will process
        startIdx = (cc-1)*chunkSize + 1;
        endIdx = min([cc*chunkSize nGood]);
        subChunk = (endIdx - startIdx) + 1;

        % Clear the par return vars
        par_p = cell(subChunk,1); par_fVal = cell(subChunk,1);
        par_Y = cell(subChunk,1); par_yFit = cell(subChunk,1);

        % Open the parpool (and do so silently within an evalc)
        %evalc('gcp()');

        % Loop over the valid voxels
        for gg = 1:subChunk

            tic;

            % The idx of the goodIdx set that we will process
            thisIdx = (gg-1)+startIdx;

            % Loop through the stimuli and frequencies and
            % obtain the data
            rawVals = [];
            adjustedVals = [];
            for dd = 1:nStims
                for ff = 1:length(allFreqs)
                    subString = sprintf(['f%dHz_' stimulusDirections{dd}],allFreqs(ff));
                    idx = find(contains(stimLabels,subString));
                    rawVals(dd,ff,:) = results.params(goodIdx(thisIdx),idx);
                end

                % Adjust the values for the zero frequency
                for ff = 2:length(allFreqs)
                    adjustedVals(dd,ff-1,:) = rawVals(dd,ff,:)-rawVals(dd,1,:);
                end
            end

            % Set our data and weights for this vertex
            W = 1./std(adjustedVals,0,3,"omitmissing");
            Y = mean(adjustedVals,3,"omitmissing");

            % Fitting the L-M stimuli influence both the midget and parasol
            % populations; we give these weights a bit of extra oomph here
            % to encourage a good match to the data.
            W(1,:,:) = W(1,:,:)*2;

            % Interpolate the data and weights to the intermediate temporal
            % frequencies. This is to avoid weird, wigly over-fits to the
            % data
            Yinterp = []; Winterp = [];
            for whichStim = 1:nStims
                Yinterp(whichStim,:) = interp1(1:nFreqs,squeeze(Y(whichStim,:)),1:0.5:nFreqs);
                Winterp(whichStim,:) = interp1(1:nFreqs,squeeze(W(whichStim,:)),1:0.5:nFreqs);
            end

            % Get the eccentricity of this vertex
            eccDeg = eccenMap(goodIdx(thisIdx));

            % The bistratified model currently fails for eccentricities
            % below 0.5 degrees. We keep the eccentricity in range here
            eccDegThresh = max([eccDeg 0.501]);

            % Get the rfRetinal symbolic equations for this eccDeg. This
            % allows us to pass these into the BADs search and not
            % re-compute on each iteration, saving time.
            [~, ~, rfRetinal] = returnTTFAtEcc(p0,stimulusDirections,eccDegThresh,interpFreqs);

            % Define the function to return the TTF matrix
            myResponseMatrix = @(p) returnTTFAtEcc(p,stimulusDirections,eccDegThresh,interpFreqs,rfRetinal,[]);

            % Define the objective
            myObj = @(p) objFunc(p,Winterp,Yinterp,myResponseMatrix);

            % Turn off the warning
            warning('off','bads:pbUnspecified');

            % BADS it
            [p, fVal] = bads(myObj,p0,lb,ub,[],[],[],optionsBADS);

            % If the fVal is greater than the re-fit threshold, refit from
            % a different p0
            refit = false;
            fValInitial = [];
            pInitial = [];
            if fVal > refitThresh
                fValInitial = fVal;
                pInitial = p;
                [p, fVal] = bads(myObj,p1,lb,ub,[],[],[],optionsBADS);
                if fVal > fValInitial
                    fVal = fValInitial;
                    p = pInitial;
                end
                refit = true;
            end

            % Get the fit at the plotting frequencies
            myResponseMatrix = @(p) returnTTFAtEcc(p,stimulusDirections,eccDegThresh,freqsForPlotting,rfRetinal,[]);
            yFit = myResponseMatrix(p);

            % Store this result in a structure
            par_p{gg} = p;
            par_fVal{gg} = fVal;
            par_Y{gg} = Y;
            par_yFit{gg} = yFit;
            par_eccDeg{gg} = eccDeg;
            par_polarAngle{gg} = polarMap(goodIdx(thisIdx));
            par_area{gg} = vArea(goodIdx(thisIdx));

            % Send the news to the console
            elapsedMins = toc()/60;
            if refit
                fprintf('vertex %d / %d complete in %2.1f mins; fVal = %2.2f [%2.2f before refit] \n',thisIdx,nGood,elapsedMins,fVal,fValInitial);
            else
                fprintf('vertex %d / %d complete in %2.1f mins; fVal = %2.2f \n',thisIdx,nGood,elapsedMins,fVal);
            end

        end

        % Place the par loop results into a results file
        fitResults.p(goodIdx(startIdx:endIdx),:) = cell2mat(par_p(1:subChunk));
        fitResults.fVal(goodIdx(startIdx:endIdx)) = cell2mat(par_fVal(1:subChunk));
        fitResults.Y(goodIdx(startIdx:endIdx)) = par_Y(1:subChunk);
        fitResults.yFit(goodIdx(startIdx:endIdx)) = par_yFit(1:subChunk);
        fitResults.eccDeg(goodIdx(startIdx:endIdx)) = cell2mat(par_eccDeg(1:subChunk));
        fitResults.polarAngle(goodIdx(startIdx:endIdx)) = cell2mat(par_polarAngle(1:subChunk));
        fitResults.area(goodIdx(startIdx:endIdx)) = cell2mat(par_area(1:subChunk));

        % Save the results file for this subject
        filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
        save(filePath,'fitResults')

        % Shut down the parpool to prevent pseudo memory leaks that arise
        % from the persistence of symbolic toolbox elements in the workers
        evalc('delete(gcp(''nocreate''))');

    end

end

% Restore the warnstate
warning(warnstate);


%% LOCAL FUNC

function fVal = objFunc(p,Winterp,Yinterp,responseFunc)

% Use the passed function to obtain the TTFs
responseMatrix = responseFunc(p);

% The weighted fit error
fVal = norm(Winterp(:).* (Yinterp(:) - responseMatrix(:)));

% Add a penalty if the TTF is rising at the lowest frequency. This is
% invariably associated with a bad fit; we don't see vertices with this
% behavior in practice.
fVal = fVal + sum(diff(responseMatrix(:,1:2),1,2)<0)*1e6;


end
