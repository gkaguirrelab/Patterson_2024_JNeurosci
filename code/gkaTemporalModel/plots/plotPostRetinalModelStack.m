% scriptCreatePlots

% Housekeeping
clear
%close all



% Params that control the plot appearance
spacing = 5000;
stimOrder = [2 3 1];

plotColor={'k','k','r','b'};
figure
figuresize(500,300,'pt');


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
            stimulusContrastScale = returnStimulusContrastScale(cellClassIndices{cc},stimulusDirections{ss});

            % Obtain this temporal RF
            rfRGC = returnPostRetinalRF(cellClassIndices{cc},stimulusDirections{ss},rgcTemporalModel,eccDeg,stimulusContrastScale);
            yVals = abs(double(subs(rfRGC,myFreqs)));

            subplot(1,4,subIdx(ss));
            offset = 1+(nEccBands*spacing) - (ee)*spacing;

            % Show the data itself
            semilogx(myFreqs,yVals+offset,...
                lineSpec{ss},'Color',plotColor{subIdx(ss)});
            hold on

            % Add a refline
            semilogx([1 1],[offset offset+2000],'-k');
            semilogx([1 100],[offset offset],':','Color',[0.5 0.5 0.5]);

            % Add a text label for the eccentricitiy
            text(80,offset+spacing/2,sprintf('%2.0f°',eccDeg));


        end
    end
end

% Clean up
for ss=1:4
    subplot(1,4,ss)
    semilogx([16 16],[-25000 50000],'--k')
    xlim([0.5 150])
    ylim([-25000 50000])
    a=gca;
    a.XTick = [1,100];
    a.XTickLabel = {'1','100'};
    a.XTickLabelRotation = 0;
    a.XMinorTick = 'off';
    a.YAxis.Visible = 'off';
    box off
end


% Save the plot
plotName = 'postRetinalStackedModelFits.pdf';
saveas(gcf,fullfile(savePath,plotName));

