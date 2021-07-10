%% makeTTFPlots
%
%

analysisIDs = { {'60d79805bdc454c70ce074ba','60d797efb57eec5d8896ca85','60d797d64137d32dd02b6dc2'} , ...
    {'60d7983aee644e04d2e077f2', '60d7982a2f9897247256ae78', '60d798187f0b6d3546dd8edf'} };

% Create a flywheel object
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download, unzip, and load the retino maps
retinoMapID = '5dc88aaee74aa3005e169380';
retinoFileName = 'TOME_3021_cifti_maps.zip';
[~, userName] = system('whoami');
userName = strip(userName);        
scratchSaveDir = ['/Users/' userName '/Desktop/tempFiles'];
saveDir = fullfile(scratchSaveDir,'v0','output');
tmpPath = fullfile(saveDir,retinoFileName);
fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
command = ['unzip -q -n ' tmpPath ' -d ' saveDir];
system(command);

% Load the retino maps
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_varea.dtseries.nii'));
vArea = cifti_read(tmpPath);
vArea = vArea.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_eccen.dtseries.nii'));
eccenMap = cifti_read(tmpPath);
eccenMap = eccenMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_angle.dtseries.nii'));
polarMap = cifti_read(tmpPath);
polarMap = polarMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_sigma.dtseries.nii'));
sigmaMap = cifti_read(tmpPath);
sigmaMap = sigmaMap.cdata;

% Create a "subcortical" map
subcorticalMap = zeros(size(vArea));
subcorticalMap(1:26298)=1;

% Eccentricity divider
eccenDivide = [0,2.5,5,10,20,40];
lineColors = {[1 0 0],[0.75 0.25 0.25],[0.5 0.5,0.5],[0.25 0.25 0.75],[0 0 1]};

eccenDivide = [0,15,30,60];
eccenDivide = [0,60];
lineColors = {[1 0 0],[0.5 0.5,0.5],[0 0 1]};


% Save location for the maps
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'LF','L-M','S'};

% Loop over subjects
for ss = 1: length(subjectNames)
    
    % Set up the paths for this subject
    fileStem = [subjectNames{ss} '_agtcOL_'];
    resultsSaveDir = ['/Users/aguirre/Desktop/' subjectNames{ss}];
    
    % Set up the figure
    figHandle = figure();
    orient(figHandle,'landscape')
    
    % Loop over analysis IDs
    for aa = 1:length(analysisIDs{ss})
        
        subplot(1,length(analysisIDs{ss}),aa);
        
        % Load the TTF model fits
        outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} '_ResultsFit.mat']);
        load(outFile,'resultsFit');
        
        vertexCount = [];
        
        for ee = 1:length(eccenDivide)-1
            
            % Get the vertices in this eccentricity band
            goodIdx = logical( (resultsFit.R2 > 0.15) .* (resultsFit.fitR2 > 0.8) .* (eccenMap > eccenDivide(ee)) .* (eccenMap < eccenDivide(ee+1)) .* (or(vArea==1,vArea==1)) );
%            goodIdx = logical( (resultsFit.R2 > 0.15) .* (resultsFit.fitR2 > 0.8) .* (subcorticalMap==1) );
            scaleMax = max(nanmean(resultsFit.fitStuff.myFit(goodIdx,:)));
            semilogx(resultsFit.fitSupport.freqs,nanmean(resultsFit.fitStuff.yVals(goodIdx,:))./scaleMax,'*','Color',lineColors{ee});
            hold on
            semilogx(resultsFit.fitSupport.freqsFit,nanmean(resultsFit.fitStuff.myFit(goodIdx,:))./scaleMax,'-','Color',lineColors{ee})
            vertexCount(ee) = sum(goodIdx);
        end
        
        % Clean up
        title(sprintf([shortNames{ss} '-' analysisLabels{aa} ' v = [%d,%d,%d]'],vertexCount))
        ylim([-.1 1.1]);
        ylabel('BOLD response [% change]')
        xlabel('Frequency [hz]');
    end
    
    % Save the figure
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_TTF.pdf']);
    print(figHandle,outFile,'-dpdf','-bestfit');
end
