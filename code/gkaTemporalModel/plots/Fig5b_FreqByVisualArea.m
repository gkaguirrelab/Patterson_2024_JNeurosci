

%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;
nBoots = 200;

% These variables define the subject names, stimulus directions. The
% Flywheel analysis IDs are listed for completeness, but not used here.
% Other software downloads the files from Flywheel.
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(studiedFreqs);
interpFreqs = logspace(log10(1),log10(100),501);
nAcqs = 12;

% Define some ROI sets
roiSet = {'LGN','V1','V2/V3','hV4','MT'};
nROIs = length(roiSet);

%% Download Mt Sinai results
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies. The response for each acquisition is retained to support
% subsequent boot-strap resampling of the data.

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('mriSinaiAnalysis'),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the LGN ROI
tmpPath = fullfile(localDataDir,'retinoFiles','LGN_bilateral.dtseries.nii');
LGNROI = cifti_read(tmpPath); LGNROI = LGNROI.cdata;

% Load the MT ROI
tmpPath = fullfile(localDataDir,'retinoFiles','MT.dtseries.nii');
MTROI = cifti_read(tmpPath); MTROI = MTROI.cdata;

% Prepare the figure
figHandle = figure('Renderer','painters');
figuresize(600,400,'pt');
tiledlayout(1,2,'TileSpacing','tight','Padding','tight')

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};


%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Loop over bootstraps
    parfor bb = 1:nBoots

        % Get a sampling (with replacement) of the 12 acquisitions
        bootIdx = datasample(1:nAcqs,nAcqs);

        % Define some variables for parpool happiness
        nGood = []; peakFreq = []; peakAmp = []; goodIdx = [];

        for rr = 1:length(roiSet)

            switch roiSet{rr}
                case 'LGN'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (LGNROI == 1)));
                case 'V1'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 1)));
                case 'V2/V3'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea >= 2) .* (vArea <= 3) ));
                case 'hV4'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea == 4) ));
                case 'MT'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (MTROI == 1)));
            end
            nGood(rr) = length(goodIdx);

            % Loop over the stimuli
            for whichStim = 1:nStims

                % Loop through the stimuli and frequencies and obtain the data
                rawVals = [];
                adjustedVals = [];

                for ff = 1:length(allFreqs)
                    subString = sprintf(['f%dHz_' stimulusDirections{whichStim}],allFreqs(ff));
                    idx = find(contains(stimLabels,subString));
                    rawVals(ff,:,:) = results.params(goodIdx,idx);
                end

                % Adjust the values for the zero frequency
                for ff = 2:length(allFreqs)
                    adjustedVals(ff-1,:,:) = rawVals(ff,:,:)-rawVals(1,:,:);
                end

                % Take the mean across voxels
                adjustedVals = squeeze(mean(adjustedVals,2));

                % Get the mean across (bootstrap resampled) acquisitions
                W = 1./std(adjustedVals(:,bootIdx),0,2)';
                Y = mean(adjustedVals(:,bootIdx),2)';

                [p,~,~,yFitInterp] = fitWatsonModel(Y,W,studiedFreqs);

                % Determine the peak frequency in the log domain
                peakFreq(whichStim,rr) = log10(interpFreqs(yFitInterp==max(yFitInterp)));

                % Save the peak amplitude, which is given by the first
                % param value
                peakAmp(whichStim,rr) = p(1);

            end % stimuli

        end % ROIs

        % Store the bootstrap result in a par cell variable
        par_peakFreq{bb} = peakFreq;
        par_peakAmp{bb} = peakAmp;

    end % Bootstraps

    % Reshape the par cell data into matrices and get mean and SEM
    peakFreq =  reshape(cell2mat(par_peakFreq),nStims,nROIs,nBoots);
    peakFreqIQR(ss,:,:) = 10.^iqr(peakFreq,3);
    peakFreqMedian(ss,:,:) = 10.^median(peakFreq,3);

    peakAmp =  reshape(cell2mat(par_peakAmp),nStims,nROIs,nBoots);
    peakAmpIQR(ss,:,:) = iqr(peakAmp,3);
    peakAmpMedian(ss,:,:) = median(peakAmp,3);

    % Plot the results
    nexttile()

    for whichStim  = 1:nStims
        vec = squeeze(peakFreqMedian(ss,whichStim,:))';
        veciqr = squeeze(peakFreqIQR(ss,whichStim,:))';
        x = (1:nROIs) + (stimOrder(whichStim)-2)/4;
        for rr=1:length(x)
            plot([x(rr) x(rr)],[vec(rr)-veciqr(rr)/2, vec(rr)+veciqr(rr)/2],'-','Color',plotColor{stimOrder(whichStim)});
        hold on
        end
        plot(x,vec,'o','Color',plotColor{stimOrder(whichStim)});
        plot(x,vec,':','Color',plotColor{stimOrder(whichStim)});
    end

            a = gca();
    a.XTick = 1:nROIs;
    a.XTickLabel = roiSet;

    title(subjects{ss});
    drawnow

end


% Save the plot
plotNamesPDF = 'Fig5b_peakFreqByArea.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

