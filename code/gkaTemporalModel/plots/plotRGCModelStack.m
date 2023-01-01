% scriptCreatePlots

% Housekeeping
clear
%close all



% Params that control the plot appearance
spacing = 200;
stimOrder = [2 3 1];

plotColor={[0.75 0.75 0.75],[0.75 0.75 0.75],[0.75 0.5 0.5],[0.5 0.5 0.75]};
lineColor={'k','k',[.5 0.25 0.25],[0.25 0.25 0.5]};

figure('Renderer','painters');
figuresize(800,600,'pt');
set(gcf, 'Color', 'None');
tiledlayout(1,4,'TileSpacing','none','Padding','tight')


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
myFreqs = logspace(log10(1),log10(100),101);
myEccs = [1,2,4,8,16,32,64];

xOffsets = (length(myEccs):-1:1)/4;

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
                subIdx=[3,2];
                lineSpec = {'-','-'};
            case 'parasol'
                stimulusDirections = {'LMS'};
                subIdx=1;
                lineSpec = {'-'};
            case 'bistratified'
                stimulusDirections = {'S'};
                subIdx=4;
                lineSpec = {'-'};
        end

        % Loop over the stimulus directions
        for ss = 1:length(stimulusDirections)

            % Obtain the chromatic weights
            [chromaticCenterWeight,chromaticSurroundWeight] = ...
                returnRGCChromaticWeights(cellClassIndices{cc},stimulusDirections{ss},eccDeg,LMRatio);

            % Obtain this temporal RF
            rfRGC = returnRGCRF(pRGCBlock,cfCone,coneDelay,chromaticCenterWeight,chromaticSurroundWeight);
            yVals = abs(double(subs(rfRGC,myFreqs)));

            % Get the subplot and define the offset
            nexttile(subIdx(ss));
            set(gca, 'XScale', 'log')
            yOffset = 1+(nEccBands*spacing) - (ee)*spacing;

            % Create a patch for the response
            X = [myFreqs fliplr(myFreqs)];
            Y = repmat(ee*spacing,size(X));
            Z = [yVals, zeros(size(yVals))];
            p = fill3(X,Y,Z,plotColor{subIdx(ss)});
            set(p,'edgecolor','none','facealpha',1);
            hold on

            % Add a line at the edge of the color patch
            plot3(myFreqs,repmat(ee*spacing,size(myFreqs)),yVals,...
                lineSpec{ss},'Color',lineColor{subIdx(ss)},'LineWidth',1.5);

            % Add a text label for the eccentricitiy
            text(150,ee*spacing,0,sprintf('%2.0f°',eccDeg),"HorizontalAlignment","left");

        end
    end
end

% Clean up
for ss=1:4
    view(0,70);
    nexttile(ss)
    xlim([0.5 150])
    zlim([0 10])
    a=gca;
    a.XTick = [1,100];
    a.XTickLabel = {'1','100'};
    a.XTickLabelRotation = 0;
    a.XMinorTick = 'off';
    a.YAxis.Visible = 'off';
    box off
    if ss>1
        a.XAxis.Visible = 'off';
        a.ZAxis.Visible = 'off';
    end

end


% Save the plot
% Save just the vector elements
for ss=1:4
    nexttile(ss);
    child_handles = allchild(gca);
    for cc=1:length(child_handles)
        if isgraphics(child_handles(cc),'patch')
            set(child_handles(cc), 'Visible','off')
        else
            set(child_handles(cc), 'Visible','on')
        end
    end
end
saveas(gcf,fullfile(savePath,'rgcStackedModelFits.pdf'));

for ss=1:4
    nexttile(ss);

    child_handles = allchild(gca);
    for cc=1:length(child_handles)
        if ~isgraphics(child_handles(cc),'patch')
            set(child_handles(cc), 'Visible','off')
        else
            set(child_handles(cc), 'Visible','on')
        end
        a=gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        a.ZAxis.Visible = 'off';
    end
end
print(gcf,fullfile(savePath,'rgcStackedModelFits.png'),'-dpng','-r600');


