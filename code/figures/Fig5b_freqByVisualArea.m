

%% Housekeeping
clear
close all
rng(1); % Fix the random number generator
verbose = false;

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

nBoots = 100;

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
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
roiSet = {'V1','V2/V3','hV4','V3a/b'};
nROIs = length(roiSet);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% Prepare the figure
figHandle = figure('Renderer','painters');
figuresize(600,600,'pt');

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
subMarkers = {'^','square','o'};
subLines = {'-','--',':'};
dirPlotShift = 0.1;


%% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Get the "Benson" visual areas for this subject
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    vAreas = cifti_read(tmpPath); vAreas = vAreas.cdata;

    % Get the "Wang" visual areas for this subject
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_wang.dscalar.nii']);
    vWang = cifti_read(tmpPath); vWang = vWang.cdata;

    % Get the area MT ROI
    tmpPath = fullfile(localDataDir,'MT.dtseries.nii');
    mtROI = cifti_read(tmpPath); mtROI = mtROI.cdata;

    % Loop over bootstraps
    parfor bb = 1:nBoots

        % Get a sampling (with replacement) of the 12 acquisitions
        bootIdx = datasample(1:nAcqs,nAcqs);

        % Define some variables for parpool happiness
        nGood = []; peakFreq = []; peakAmp = []; goodIdx = [];

        % Loop over the ROIs
        for rr = 1:length(roiSet)

            switch roiSet{rr}
                case 'V1'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas == 1)));
                case 'V2/V3'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas >= 2) .* (vAreas <= 3) ));
                case 'hV4'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas == 4) ));
                case 'V3a/b'
                    goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas >= 11) .* (vAreas <= 12) ));
                otherwise
                    error('I do not know that region')
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

                % Determine the peak frequency in the log domain. There's a
                % little business here to handle the edge case of more than
                % one frequency being equally maximal.
                peakFreq(whichStim,rr) = mean(log10(interpFreqs(yFitInterp==max(yFitInterp))));

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

    for whichStim  = 1:nStims
        vec = squeeze(peakFreqMedian(ss,whichStim,:))';
        veciqr = squeeze(peakFreqIQR(ss,whichStim,:))';
        x = (1:nROIs);
        x = x + (ss-2)*dirPlotShift;
        plot(x,vec,subMarkers{ss},...
            'MarkerEdgeColor','none','MarkerFaceColor',plotColor{stimOrder(whichStim)},...
            'MarkerSize',10)
        hold on
        if whichStim == 3
            plot(x,vec,':','Color',plotColor{stimOrder(whichStim)},'LineWidth',1);
        end

        % Save the vector for the across-subject plot
        bigVec(ss,whichStim,:)=vec;
    end

    a = gca();
    a.XTick = 1:nROIs;
    a.XTickLabel = roiSet;

    ylim([0 50])
    drawnow

end

x = 1:nROIs;
for whichStim  = 1:nStims
    y = median(squeeze(bigVec(:,whichStim,:)));
    pH(whichStim) = plot(x,y,'-','Color',plotColor{stimOrder(whichStim)},'LineWidth',2);
    plot(x,y,'o',...
        'MarkerFaceColor','none','MarkerEdgeColor',plotColor{stimOrder(whichStim)},...
        'MarkerSize',20,'LineWidth',2)

end
legend(pH,subjects)

% Save the plot
plotNamesPDF = 'Fig5b_peakFreqByArea.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));

