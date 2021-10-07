%% demoExtractTTFData
%
% This script demonstrates how to download the "results" files Flywheel and
% extract BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
shortNames = {'gka','asb'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
nFreqs = length(freqs);

% Frequency components for model fitting
deltaF10 = min(diff(log10(freqs(2:end))));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);

% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the subcortical ROIs
projectID = '5ca7803af546b60029ef118e';
subCorticalROIsFullNames = {'LGN_bilateral.dtseries.nii','thalamus_bilateral.dtseries.nii','midbrain_bilateral.dtseries.nii'};
subCorticalROIsLabels = {'LGN','thalamus','midbrain'};
for rr = 1:length(subCorticalROIsFullNames)
    tmpPath = fullfile(localSaveDir,'retinoFiles',subCorticalROIsFullNames{rr});    
    fw.downloadFileFromProject(projectID,subCorticalROIsFullNames{rr},tmpPath);
    tmpRegion = cifti_read(tmpPath); tmpRegion = tmpRegion.cdata;
    str = [subCorticalROIsLabels{rr} 'ROI = tmpRegion;'];
    eval(str);
end

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.1;

% This is the visual area and eccentricity range to grab. 
areaLabel = 'LGN';
eccenRange = [0 90];

% Define some components for model fitting
p0 = [2, 0.9, 0.015];
lb = [0 0 0];
ub = [10 1 0.025];

options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% Create a figure
figure;

% Create a data variable to hold the results. This will be a 2 x 3 x 6
% (subjects x directons x frequencies) cell array for the selected subject. Each cell
% will have the 12 measurements
data = cell(2,3,6);

% Loop through the directions
for ss = 1:2
    for dd = 1:3
        
        % Download the results file for this subject        
        fileName = [subjectNames{ss} '_mtSinai_results.mat'];
        filePath = fullfile(localSaveDir,'resultsFiles',fileName);
        fw.downloadOutputFromAnalysis(analysisIDs{ss},fileName,filePath);
        load(filePath,'results')
        
        % Grab the stimLabels
        stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};
        
        % Find the vertices that we wish to analyze
        switch areaLabel
            case 'V1'
                areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'V23'
                areaIdx = (vArea>=2) .* (vArea<=3) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'hV4'
                areaIdx = (vArea>=4) .* (vArea<=5) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'MT'
                areaIdx = (vArea>=8) .* (vArea<=9) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'LGN'
                areaIdx = LGNROI;
            case 'midbrain'
                areaIdx = midbrainROI;
            case 'thalamus'
                areaIdx = thalamusROI;
        end
        
        goodIdx = logical( (results.R2 > r2Thresh) .* areaIdx );
        
        % Loop through the frequencies and obtain the set of values
        vals = cell(1,nFreqs);
        for ff = 1:nFreqs
            subString = sprintf(['f%dHz_' directions{dd}],freqs(ff));
            idx = find(contains(stimLabels,subString));
            vals{ff} = mean(results.params(goodIdx,idx),'omitnan');
        end
        
        % Prepare to plot into this subplot
        subplot(2,3,dd+(ss-1)*3);
        
        % Adjust the values for the zero frequency and plot
        for ff = 2:nFreqs
            data{ss,dd,ff-1} = vals{ff}-vals{1};
            semilogx(zeros(1,length(data{ss,dd,ff-1}))+freqs(ff),data{ss,dd,ff-1},'.','Color',[0.5 0.5 0.5]);
            hold on
        end
        
        % Obtain the mean across frequencies and add to the plot
        meanVals = cellfun(@(x) mean(x,'omitnan'),squeeze(data(ss,dd,:)));
        semilogx(freqs(2:end),meanVals,'ob','MarkerSize',10,'MarkerFaceColor','b')
        
        % Fit the Watson model
        myObj = @(p) norm(meanVals' - watsonTTF(p,freqs(2:end)));
        p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
        myFit = watsonTTF(p,freqsFit);
        
        % Add the fitted TTF
        semilogx(freqsFit,myFit,'-r');
        
        % Clean up the plot
        title([shortNames{ss} ' ' directions{dd} sprintf('n=%d',sum(goodIdx))])
        ylabel('BOLD % change');
        xlabel('frequency [Hz]');
        semilogx([1 64],[0 0],':k','LineWidth',1)
        ylim([-2 8]);
        xlim([1 128])
        set(gca,'TickDir','out');
        box off
        
    end
end
