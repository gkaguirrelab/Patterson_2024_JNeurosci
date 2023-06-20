% plotRGCModelComponents

% Where to save figures
savePath = fullfile('~','Desktop','VSS 2023');

% Load the RGC temporal model
rgcTemporalModel = fitRGCFResponse();

% Pick an eccentricity, cell, and stimulus direction
eccDeg = 20.5;
cellClass = 'midget';
stimulusDirection = 'LMS';
stimulusContrastScale = returnStimulusContrastScale(cellClass,stimulusDirection);
[rfPostRetinal, rfRGC, rfBipolar, rfCone] = returnRetinalRF(cellClass,stimulusDirection,rgcTemporalModel,eccDeg,stimulusContrastScale);

% Plot the cell IRFs
figHandle = figure();
figuresize(200, 400,'pt');
plotRF(rfCone,figHandle,'-g');
plotRF(rfBipolar,figHandle,'-m');
plotRF(rfRGC,figHandle,'-c');

% Save the figure
plotName = [cellClass '_' stimulusDirection '_' num2str(eccDeg,2) '_CellIRFs.pdf' ];
saveas(figHandle,fullfile(savePath,plotName));

% Now plot the filters
cellIndex = strcmp(cellClass,{'midget','parasol','bistratified'});
pRGCBlock = nan(1,7);
for ii = 1:7
    pRGCBlock(ii) = rgcTemporalModel.pFitByEccen{ii,cellIndex}(eccDeg);
end

% Extract the parameters
g = pRGCBlock(1); k = pRGCBlock(2);
cfInhibit = pRGCBlock(3); cf2ndStage = pRGCBlock(4); Q = pRGCBlock(5);
surroundWeight = pRGCBlock(6); surroundDelay = pRGCBlock(7);

figHandle = figure();
figuresize(200, 400,'pt');

filter = stageInhibit(cfInhibit,k);
plotRF(filter,figHandle,'-g');

filter = stageSecondOrderLP(cf2ndStage,Q);
plotRF(filter,figHandle,'-k');

plotName = [cellClass '_' stimulusDirection '_' num2str(eccDeg,2) '_FilterIRFs.pdf' ];
saveas(figHandle,fullfile(savePath,plotName));
