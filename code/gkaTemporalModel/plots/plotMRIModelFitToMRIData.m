% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Extract some info from the stored model structure
pRGC = rgcTemporalModel.p;
cfCone = rgcTemporalModel.cfCone;
coneDelay = rgcTemporalModel.coneDelay;
LMRatio = rgcTemporalModel.LMRatio;
pFitByEccen = rgcTemporalModel.pFitByEccen;
blockParamNames = rgcTemporalModel.meta.blockParamNames;
eccFields = rgcTemporalModel.meta.eccFields;
eccBins = rgcTemporalModel.meta.eccBins;
lbBlock = rgcTemporalModel.meta.lbBlock;
ubBlock = rgcTemporalModel.meta.ubBlock;
eccDegs = rgcTemporalModel.meta.eccDegs;

nEccBands = length(eccFields);
nBlockParams = size(pRGC,1);



%% plot the model fits to the MRI data

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','mriTemporalModel.mat');
load(loadPath,'mriTemporalModel');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriTemporalModel.meta.studiedFreqs;
eccDegVals = mriTemporalModel.meta.eccDegVals;
subjects = mriTemporalModel.meta.subjects;
stimulusDirection = mriTemporalModel.meta.stimulusDirection;
plotColor = mriTemporalModel.meta.plotColor;
nFixed = mriTemporalModel.meta.nFixed;
nEcc = length(eccDegVals);
freqsForPlotting = logspace(0,2,50);

for whichStim = 1:length(stimulusDirection)
    for whichSub = 1:length(subjects)

        figure
        pMRI = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).p;
        v1Y = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Y;
        v1Y = reshape(v1Y,length(studiedFreqs),1,nEcc);
        v1Eccentricity = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.v1Eccentricity;
        lgnFreqX = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnFreqX;
        lgnY = mriTemporalModel.(stimulusDirection{whichStim}).(subjects{whichSub}).data.lgnY;

        % Loop over eccentricities
        for ee=1:nEcc
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1Y(:,1,ee)),['o' plotColor{whichStim}]);
            hold on
            pMRIBlock = [pMRI(1:5) pMRI(nFixed+ee) pMRI(nFixed+nEcc+ee)];
            yFit = returnV1TTFForEcc(stimulusDirection{whichStim},rgcTemporalModel,eccDegVals(ee),pMRIBlock,freqsForPlotting);
            semilogx(freqsForPlotting,yFit,['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimulusDirection{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
            ylim([-1 7]);
        end

        % Add the LGN response
        lgnTTFFit = returnlgnTTF(stimulusDirection{whichStim},rgcTemporalModel,pMRI,freqsForPlotting,v1Eccentricity);
        subplot(2,4,8)
        semilogx(lgnFreqX,lgnY,['o' plotColor{whichStim}]);
        hold on
        semilogx(freqsForPlotting,lgnTTFFit,['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimulusDirection{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

        % Plot the surround suppression index vs. eccentricity
        subplot(2,4,4)
        plot(log10(eccDegVals),pMRI(:,nFixed+1:nFixed+nEcc),['*' plotColor{whichStim}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Save the plot
        plotName = [stimulusDirection{whichStim} '_' subjects{whichSub} '_ModelFit.pdf' ];
        saveas(gcf,fullfile('~/Desktop',plotName));
    end
end
