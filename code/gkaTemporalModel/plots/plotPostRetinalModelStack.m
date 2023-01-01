% scriptCreatePlots

% Housekeeping
clear
%close all



% Params that control the plot appearance
spacing = 200;
stimOrder = [2 3 1];

plotColor={'k',[0.75 0.5 0.5],[0.5 0.5 0.75]};
lineColor={'none',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = [0,1,1];
    figure('Renderer','painters');
figuresize(600,300,'pt');


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

myFreqs = logspace(log10(1),log10(100),101);
myEccs = [1,2,4,8,16,32,64];

xOffsets = (length(myEccs):-1:1)/4;

for ee = 1:length(myEccs)

    eccDeg = myEccs(ee);

    totalAchrom = zeros(size(myFreqs));
    achromCounter = 0;

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
                subIdx=[2,1];
                lineSpec = {'-',':'};
            case 'parasol'
                stimulusDirections = {'LMS'};
                subIdx=1;
                lineSpec = {'-'};
            case 'bistratified'
                stimulusDirections = {'S'};
                subIdx=3;
                lineSpec = {'-'};
        end

        % Loop over the stimulus directions
        for ss = 1:length(stimulusDirections)

            % Obtain the stimulus contrast weights
            stimulusContrastScale = returnStimulusContrastScale(cellClassIndices{cc},stimulusDirections{ss});

            % Obtain this temporal RF
            rfRGC = returnPostRetinalRF(cellClassIndices{cc},stimulusDirections{ss},rgcTemporalModel,eccDeg,stimulusContrastScale);

            % Get the yVals and scale them to simplify plotting
            yVals = abs(double(subs(rfRGC,myFreqs)));
            yVals = yVals ./ 1e4;

            % Increase the bistratifieds so they can be seen
            switch stimulusDirections{ss}
                case 'LMS'
                    achromCounter = achromCounter+1;
                    totalAchrom = totalAchrom + yVals;
                    yVals = yVals;
                case 'S'
                    yVals = yVals * 50;
            end

            % Get the subplot and define the offset
            subplot(1,3,subIdx(ss));
            set(gca, 'XScale', 'log')

            % Create a patch for the response
            X = [myFreqs fliplr(myFreqs)];
            Y = repmat(ee*spacing,size(X));
            Z = [yVals, zeros(size(yVals))];
             p = fill3(X,Y,Z,plotColor{subIdx(ss)});
             set(p,'edgecolor','none','facealpha',faceAlpha(subIdx(ss)));
            hold on

            % Add a line at the edge of the color patch            
            plot3(myFreqs,repmat(ee*spacing,size(myFreqs)),yVals,...
                lineSpec{ss},'Color',lineColor{subIdx(ss)},'LineWidth',1.5);

            % Add a text label for the eccentricitiy
            text(150,ee*spacing,0,sprintf('%2.0f°',eccDeg),"HorizontalAlignment","left");

            % Plot the total achrom case
            if achromCounter==2
                subplot(1,3,1);
                set(gca, 'XScale', 'log')
                X = [myFreqs fliplr(myFreqs)];
                Y = repmat(ee*spacing,size(X));
                Z = [totalAchrom, zeros(size(totalAchrom))];
                p = fill3(X,Y,Z,[0.75 0.75 0.75]);
                set(p,'edgecolor','none','facealpha',1);
                hold on

                % Add a line at the edge of the color patch
                plot3(myFreqs,repmat(ee*spacing,size(myFreqs)),totalAchrom,...
                    '-','Color','k','LineWidth',1.5);

                % Add a text label for the eccentricitiy
                text(150,ee*spacing,0,sprintf('%2.0f°',eccDeg),"HorizontalAlignment","left");
            end


        end
    end
end

% Clean up
for ss=1:3
    subplot(1,3,ss)
    view(0,70);
    plot3([16 16],[0 spacing*length(myEccs)],[10 10],'--k')
    xlim([0.5 150])
    zlim([0 5])
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

