%% makeTTFPlots
%
%

r2Thresh = 0.25;
fitR2Thresh = 0.8;

analysisLabels = {'LF','L-M','S'};
analysisIDs = { {'60e9ea6dbceb4c0bc9e0767e','60e9ea50a85492ed8f96cabd','60e9ea334ef89230db2b7021'} , ...
    {'60e9eabf4ef89230db2b7027', '60e9eaa1a74445f40c56b123', '60e9ea85bd00f64426dd9301'} };

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
eccenDivide = [0,20,40,80];
eccenNames = {'0-20','20-40','40-80'};
areaNames = {'LGN','V1','V2V3','hV4+','MT+'};

% Save location for the maps
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
lineColors = {[0.5 0.5,0.5],[1 0 0],[0 0 1],[0 1 0],[0 1 1]};

% Loop over subjects
for ss = 1: length(subjectNames)
    
    % Set up the paths for this subject
    fileStem = [subjectNames{ss} '_agtcOL_'];
    resultsSaveDir = ['/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Patterson_2021_EccentricityFlicker/matlabFigures/' subjectNames{ss}];
    
    % Set up the figure
    figHandle0 = figure();
    figHandle1 = figure();
    figHandle2 = figure();
    orient(figHandle1,'landscape')
    orient(figHandle2,'landscape')
    
    acrossAnalysisVertexCount = [];
    
    % Loop over analysis IDs
    for aa = 1:length(analysisIDs{ss})
        
        % Load the TTF model fits
        outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} '_ResultsFit.mat']);
        load(outFile,'resultsFit');
        
        % All responsive voxels in V1
        set(0, 'CurrentFigure', figHandle0)
        goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (vArea==1) );
        scaleMax = max(nanmean(resultsFit.fitStuff.myFit(goodIdx,:)));
%        semilogx(resultsFit.fitSupport.freqs,nanmean(resultsFit.fitStuff.yVals(goodIdx,:))./scaleMax,'*','Color',lineColors{aa});
        semilogx(resultsFit.fitSupport.freqsFit,nanmean(resultsFit.fitStuff.myFit(goodIdx,:))./scaleMax,'-','Color',lineColors{aa})
        hold on
        acrossAnalysisVertexCount(aa) = sum(goodIdx);
        
        % Loop over eccentricities within V1
        set(0, 'CurrentFigure', figHandle1)
        subplot(1,length(analysisIDs{ss}),aa);
        vertexCount = [];
        for ee = 1:length(eccenDivide)-1
            
            % Get the vertices in this eccentricity band
            goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (eccenMap > eccenDivide(ee)) .* (eccenMap < eccenDivide(ee+1)) .* (or(vArea==1,vArea==1)) );
            scaleMax = max(nanmean(resultsFit.fitStuff.myFit(goodIdx,:)));
            %            semilogx(resultsFit.fitSupport.freqs,nanmean(resultsFit.fitStuff.yVals(goodIdx,:))./scaleMax,'*','Color',lineColors{ee});
            semilogx(resultsFit.fitSupport.freqsFit,nanmean(resultsFit.fitStuff.myFit(goodIdx,:))./scaleMax,'-','Color',lineColors{ee})
            hold on
            vertexCount(ee) = sum(goodIdx);
        end
        
        % Clean up
        title(sprintf([shortNames{ss} '-' analysisLabels{aa} ' v = [%d,%d,%d]'],vertexCount))
        ylim([-.1 1.1]);
        ylabel('BOLD response [% change]')
        xlabel('Frequency [hz]');
        legend(eccenNames);
        
        % Loop over regions
        set(0, 'CurrentFigure', figHandle2)
        subplot(1,length(analysisIDs{ss}),aa);
        vertexCount = [];
        for rr = 1:length(areaNames)
            
            switch rr
                case 1
                    goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (subcorticalMap==1) );
                case 2
                    goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (vArea==1) );
                case 3
                    goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (or(vArea==2,vArea==3)) );
                case 4
                    goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (or(vArea==4,vArea==5)) );
                case 5
                    goodIdx = logical( (resultsFit.R2 > r2Thresh) .* (resultsFit.fitR2 > fitR2Thresh) .* (or(vArea==8,vArea==9)) );
            end
            
            % Get the voxels/vertices in this region
            scaleMax = max(nanmean(resultsFit.fitStuff.myFit(goodIdx,:)));
            %            semilogx(resultsFit.fitSupport.freqs,nanmean(resultsFit.fitStuff.yVals(goodIdx,:))./scaleMax,'*','Color',lineColors{ee});
            semilogx(resultsFit.fitSupport.freqsFit,nanmean(resultsFit.fitStuff.myFit(goodIdx,:))./scaleMax,'-','Color',lineColors{rr})
            hold on
            vertexCount(rr) = sum(goodIdx);
        end
        
        % Clean up
        title(sprintf([shortNames{ss} '-' analysisLabels{aa} ' v = [%d,%d,%d,%d,%d]'],vertexCount))
        ylim([-.1 1.1]);
        ylabel('BOLD response [% change]')
        xlabel('Frequency [hz]');
        legend(areaNames);
    end
    
    % Clean up figure0
    set(0, 'CurrentFigure', figHandle0)
    title(sprintf([shortNames{ss} '- All V1; v = [%d,%d,%d]'],acrossAnalysisVertexCount))
    ylim([-.1 1.1]);
    ylabel('BOLD response [% change]')
    xlabel('Frequency [hz]');
    legend(analysisLabels);
    
    % Save the figure
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_V1DirectionTTF.pdf']);
    print(figHandle0,outFile,'-dpdf','-bestfit');
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_V1EccenTTF.pdf']);
    print(figHandle1,outFile,'-dpdf','-bestfit');
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_AreaTTF.pdf']);
    print(figHandle2,outFile,'-dpdf','-bestfit');
end


