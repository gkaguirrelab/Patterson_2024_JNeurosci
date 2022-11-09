% scriptCreatePlots

% Housekeeping
clear
close all

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults');
load(fullfile(loadPath,'mriFullResultSet.mat'),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_RGCOnly');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
nFixedParams = mriFullResultSet.meta.nFixedParams;
nFloatByEccParams = mriFullResultSet.meta.nFloatByEccParams;
nUniqueParams = mriFullResultSet.meta.nUniqueParams;
postReceptoralPaths = mriFullResultSet.meta.postReceptoralPaths;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
chromAchromIndex = [1 1 2];
nParamsPerCellBlock = nFixedParams+nEccs*2;

for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1W = mriFullResultSet.(subjects{whichSub}).v1W;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
    lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
    lgnW = mriFullResultSet.(subjects{whichSub}).lgnW;
    lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
    lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

    switch whichSub
        case 1
            pMRI = [ ...
                11.6131162643, 0.0000000000, 0.0130615807, 0.1774803603, 0.0152459931, 200.0000000000, 0.7000000000, 200.0000000000, 0.7000000000, ... % lgn
                23.4646415710, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.7988156275, 0.8945436752, 1.9444416686, 2.0674452774, 1.6989941009, 3.8656232155, ... % V1 midget.LminusM
                21.8660521507, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.7434652722, 2.5777765641, 3.7304884986, 2.5607346525, 0.8132058998, 0.6042100380, ... % V1 bistratified.S
                5.1922225952, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.8850689925, 0.5641330718, 1.1357428575, 1.6942211341, 1.5819410767, 2.2302426601, ... % V1 parasol.LMS
                11.8660521507, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 3.6268399840, 3.3808847552, 2.1807472265, 0.3002217758, 0.3003650346 ... % V1 midget.LMS
                ];
        case 2
            pMRI = [ ...
                17.6642894745, 0.0000000000, 0.0233214378, 0.1332580185, 0.0136488342, 200.0000000000, 0.7000000000, 200.0000000000, 0.7000000000, ... % lgn
                21.1115264893, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.3998590591, 0.5808780680, 1.4821296761, 2.0052346234, 4.5578821859, 17.0652279660, ... % V1 midget.LminusM
                17.2517776489, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.8967910008, 3.9502062277, 6.0779801930, 4.6071290483, 3.7399487375, 3.9068831025, ... % V1 bistratified.S
                11.9964599609, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 1.8683113510, 1.8482511685, 0.4175473847, 2.5581598382, 2.7199083066, 4.2592832309, ... % V1 parasol.LMS
                10.5072784424, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.9463639093, 3.3549307974, 1.2737846959, 0.3719261103, 0.3001181248 ... % V1 midget.LMS
                ];
        otherwise
            % Fit the model with only RGC gain
            pMRI = fitMRIResponse(pMRI,stimulusDirections,studiedEccentricites,...
                studiedFreqs,v1Y,v1W,lgnY,lgnW,false,'rgcOnly');
            foo=1;
    end


    [~,v1YFitMatrix] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,nUniqueParams,nFixedParams);
    [~,lgnYFitMatrix] = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,freqsForPlotting,rgcTemporalModel);

    for whichStim = 1:length(stimulusDirections)
        figure
        figuresize(600,300,'pt');

        % Loop over eccentricities
        for ee=1:nEccs

            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,v1Y(v1DataIndices),['o' plotColor{whichStim}]);
            hold on

            % patch error bars
            X = [studiedFreqs fliplr(studiedFreqs)];
            Y = [v1Ylow(v1DataIndices), fliplr(v1Yhigh(v1DataIndices))];
            p = patch(X,Y,plotColor{whichStim});
            set(p,'edgecolor','none','facealpha',0.1);

            lineStyle = {'-.',':'};
            for cc=1:2
                semilogx(freqsForPlotting,squeeze(v1YFitMatrix(whichStim,ee,cc,:)),[lineStyle{cc} plotColor{whichStim}]);
            end
            semilogx(freqsForPlotting,sum(squeeze(v1YFitMatrix(whichStim,ee,:,:))),['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(studiedEccentricites(ee),2) '°']);
            ylim([-1 7]);
            a=gca; a.XTick = studiedFreqs;
            a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
            a.XTickLabelRotation = 0;
        end

        % Add the LGN response
        lgnDataIndices = 1+(whichStim-1)*(nFreqs): ...
            (whichStim-1)*(nFreqs)+nFreqs;

        subplot(2,4,8)
        semilogx(studiedFreqs,lgnY(lgnDataIndices),['o' plotColor{whichStim}]);
        hold on
        % path error bars
        X = [studiedFreqs fliplr(studiedFreqs)];
        Y = [lgnYlow(lgnDataIndices), fliplr(lgnYhigh(lgnDataIndices))];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.1);

        for cc=1:2
            semilogx(freqsForPlotting,squeeze(lgnYFitMatrix(whichStim,cc,:)),[lineStyle{cc} plotColor{whichStim}]);
        end
        semilogx(freqsForPlotting,sum(squeeze(lgnYFitMatrix(whichStim,:,:))),['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimulusDirections{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-1 7]);
        a=gca; a.XTick = studiedFreqs;
        a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
        a.XTickLabelRotation = 0;

        % Show the cortical filter
%         subplot(2,4,4)
%         secondOrderFc = pMRI(5+(chromAchromIndex(whichStim)-1)*2+1);
%         secondOrderQ = pMRI(5+(chromAchromIndex(whichStim)-1)*2+2);
% 
%         syms f; rf = stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
%         myFreqs = logspace(log10(0.5),log10(100),101);
%         ttfComplex = double(subs(rf,myFreqs));
%         gainVals = abs(ttfComplex);
%         semilogx(myFreqs,gainVals,['-.' plotColor{whichStim}]);
%         xlabel('frequency [Hz]'); ylabel('gain');
%         a=gca; a.XTick = studiedFreqs;
%         a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
%         a.XTickLabelRotation = 0;

        % Save the plot
        plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_MRIModelFit.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

    end

end
