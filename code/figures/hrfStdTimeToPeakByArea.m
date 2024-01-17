

%% Housekeeping
clear
close all
rng(1); % Fix the random number generator

% Place to save figures and to find the Watson fit results
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

%% Analysis properties
% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.1;

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};
nSubs = length(subjects);
subMarkers = {'^','square','o'};

% Define some ROI sets
roiSet = {'V1','V2/V3','hV4','V3A/B'};
nROIs = length(roiSet);

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% Get the flobs basis set with 100 ms resolution
deltaT = 0.1;
flobsbasis = returnFlobsVectors(deltaT);

figHandle = figure();
figuresize(600,300,'pt');
dirPlotShift = 0.1;

%% Loop through subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Get the "Benson" visual areas for this subject
    tmpPath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_benson.dscalar.nii']);
    vAreas = cifti_read(tmpPath); vAreas = vAreas.cdata;

    % Loop over the ROIs
    for rr = 1:length(roiSet)

        switch roiSet{rr}
            case 'V1'
                goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas == 1)));
            case 'V2/V3'
                goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas >= 2) .* (vAreas <= 3) ));
            case 'hV4'
                goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas == 4) ));
            case 'V3A/B'
                goodIdx = find(logical( (results.R2 > r2Thresh) .* (vAreas >= 11) .* (vAreas <= 12) ));
            otherwise
                error('I do not know that region')
        end
        nGood(rr) = length(goodIdx);

        % Loop over the vertices and get the time to peak of the HRF
        timeToPeak = [];
        for gg = 1:length(goodIdx)
            pp = results.params(goodIdx(gg),end-2:end);
            hrf = makeFlobsHRF(pp, flobsbasis);
            [~,idx] = max(hrf);
            timeToPeak(gg) = idx*deltaT;
        end

        % Get the mean and SD of the time-to-peak for this region
        meanTTP(ss,rr) = mean(timeToPeak);
        stdTTP(ss,rr) = std(timeToPeak);

        x = rr + (ss-2)*dirPlotShift;
        plot([x x],[meanTTP(ss,rr)+stdTTP(ss,rr),meanTTP(ss,rr)-stdTTP(ss,rr)],'-','LineWidth',2,'Color',[0.5 0.5 0.5])
        hold on
        plot(x,meanTTP(ss,rr),subMarkers{ss},'MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',12,'LineWidth',2);

    end % ROIs

end

a = gca();
a.XTick = 1:nROIs;
a.XTickLabel = roiSet;
ylim([0 8])
a.YTick = 0:2:8;
box off
ylabel('HRF time-to-peak [s]')
title('mean±SD across each region');

plotNamesPDF = 'hrfTimeToPeakByROI.pdf';
saveas(figHandle,fullfile(savePath,plotNamesPDF));




%% Local functions
function hrf = makeFlobsHRF(x, flobsbasis)

% Create the HRF
hrf = flobsbasis*x';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

end