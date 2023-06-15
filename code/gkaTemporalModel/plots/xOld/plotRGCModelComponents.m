% plotRGCModelComponents

% Where to save figures
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract some info from the stored model structure
pRGC = rgcTemporalModel.p;
cfCone = rgcTemporalModel.cfCone;
coneDelay = rgcTemporalModel.coneDelay;
LMRatio = rgcTemporalModel.LMRatio;
pFitByEccen = rgcTemporalModel.pFitByEccen;
eccFields = rgcTemporalModel.meta.eccFields;
eccBins = rgcTemporalModel.meta.eccBins;
lbBlock = rgcTemporalModel.meta.lbBlock;
ubBlock = rgcTemporalModel.meta.ubBlock;
cellClassIndices = rgcTemporalModel.meta.cellClassIndices;

nEccBands = length(eccFields);
nBlockParams = size(pRGC,1);

% Plot parasol bipolar response across eccentricity
myFreqs = logspace(log10(0.5),log10(100),101);
myEccs = 1:5:91;
for ee = 1:length(myEccs)

    eccDeg = myEccs(ee);

    for cc=1:length(cellClassIndices)

        % Obtain the RGC model parameters for this eccentricity
        pRGCBlock = [];
        for ii = 1:7
            pRGCBlock(ii) = rgcTemporalModel.pFitByEccen{ii,cc}(eccDeg);
        end

        % Determine the stim directions to plot for this cell class
        stimulusDirections = {};
        switch cellClassIndices{cc}
            case 'midget'
                stimulusDirections = {'LminusM','LMS'};
            case 'parasol'
                stimulusDirections = {'LMS'};
            case 'bistratified'
                stimulusDirections = {'S'};
        end

        % Loop over the stimulus directions
        for ss = 1:length(stimulusDirections)

            % Obtain the chromatic weights
            [chromaticCenterWeight,chromaticSurroundWeight] = ...
                returnRGCChromaticWeights(cellClassIndices{cc},stimulusDirections{ss},eccDeg,LMRatio);

            % Obtain this temporal RF
            [rfRGC, rfBipolar, rfCone] = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);

            % Construct the response grid
            plotData.(cellClassIndices{cc}).(stimulusDirections{ss}).cone(ee,:) = abs(double(subs(rfCone,myFreqs)));
            plotData.(cellClassIndices{cc}).(stimulusDirections{ss}).bipolar(ee,:) = abs(double(subs(rfBipolar,myFreqs)));
            plotData.(cellClassIndices{cc}).(stimulusDirections{ss}).rgc(ee,:) = abs(double(subs(rfRGC,myFreqs)));

        end
    end
end

stageList = {'cone','bipolar','rgc'};
stageList = {'rgc'};
zLimMax = [50,50,10];
zLimMax = [10];
cellList = {'parasol','midget','bistratified'};
stimList = {{'LMS'},{'LMS','LminusM'},{'S'}};


[X,Y] = meshgrid(log10(myFreqs),myEccs);
map = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';
map = plasmaColorMap();


for xx = 1:length(stageList)
    f(xx) = figure('Renderer','painters','Position',[100 100 600 300]);
    colormap(map)
    subPlotIndex = 1;
    for cc=1:length(cellList)
        cellType = cellList{cc};
        stimulusDirections = stimList{cc};
        for ss=1:length(stimulusDirections)
            thisData = plotData.(cellList{cc}).(stimulusDirections{ss}).(stageList{xx});
             subplot(2,4,subPlotIndex);

            s = mesh(X,Y,thisData);
            s.FaceColor = 'interp';
            s.FaceAlpha = 1;
            s.EdgeColor = 'none';
            view([15 30]);
            ylim([-5 95]);
            xlim([-0.5 2.15]);
            zlim([0 zLimMax(xx)]);
            a = gca;
            a.Color = 'none';
            a.XTickLabel={'1','10','100'};
            a.YTick=[0 30 60 90];
            a.YTickLabelRotation = 0;
            a.Color = 'none';
            a.ZTickLabel = [];
            box off
            grid off
            pbaspect([1 1 1e-6])
            title([cellList{cc} ' - ' stimulusDirections{ss}])

             subplot(2,4,subPlotIndex+4);
            s = mesh(X,Y,thisData);
            s.FaceColor = 'interp';
            s.FaceAlpha = 0.5;
            s.EdgeAlpha = 0.75;
            view([15 30]);
            ylim([-5 95]);
            xlim([-0.5 2.15]);
            zlim([0 zLimMax(xx)]);
            a = gca;
            a.Color = 'none';
            a.YTick=[0 30 60 90];
            a.YTick = [];
            a.XTick = [];
            box on
            a.BoxStyle = 'full';
            a.XColor = 'b';
            a.YColor = 'b';
            a.ZColor = 'b';
            axis square
            grid off
            subPlotIndex = subPlotIndex+1;
        end
    end
%    saveName = [stageList{xx} '_' cellList{cc} '_' stimulusDirections{ss} '_RGCModelComponents.pdf'];
%    saveas(f(xx),fullfile(savePath,saveName));
end

% 
% figure
% colormap(map)
% for ii=1:length(plotItems)
%     subplot(2,3,ii);
%     contourf(X,Y,plotData.(plotItems{ii}),15,'LineWidth',0.5)
%     axis square
% end
