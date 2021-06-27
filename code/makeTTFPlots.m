%% makeTTFPlots
%
%

% Eccentricity divider
eccenDivide = [0,2.5,5,10,20,40];
lineColors = {[1 0 0],[0.75 0.25 0.25],[0.5 0.5,0.5],[0.25 0.25 0.75],[0 0 1]};

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
        
        
        for ee = 1:length(eccenDivide)-1
            
            % Get the vertices in this eccentricity band
            goodIdx = logical( (resultsFit.R2 > 0.15) .* (resultsFit.fitR2 > 0.8) .* (eccenMap > eccenDivide(ee)) .* (eccenMap < eccenDivide(ee+1)) .* (vArea==1) );
            scaleMax = max(nanmean(resultsFit.fitStuff.myFit(goodIdx,:)));
            %        semilogx(resultsFit.fitSupport.freqs,nanmean(resultsFit.fitStuff.yVals(goodIdx,:))./scaleMax,'*','Color',lineColors{ee});
            semilogx(resultsFit.fitSupport.freqsFit,nanmean(resultsFit.fitStuff.myFit(goodIdx,:))./scaleMax,'-','Color',lineColors{ee})
            hold on
            vertexCount(ee) = sum(goodIdx);
        end
        
        % Clean up
        title(sprintf([shortNames{ss} '-' analysisLabels{aa} ' v = [%d,%d,%d,%d,%d]'],vertexCount))
        ylim([-.1 1.1]);
        ylabel('BOLD response [% change]')
        xlabel('Frequency [hz]');
    end
    
    % Save the figure
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_TTF.pdf']);
    print(figHandle,outFile,'-dpdf','-bestfit');
end
