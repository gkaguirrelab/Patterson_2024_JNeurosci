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

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_NoLateFilter');

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
16.4008557796, 0.1001822233, 0.0093419790, 0.1701574039, 0.0164101470, 200.0000000000, 0.7000000000, 200.0000000000, 0.7000000000, ... % lgn 
30.8288097382, 0.9743269920, 0.8396718740, 0.9641625643, 0.9985010624, 0.7597019911, 0.8089851141, 1.2719957935, 1.3161493631, 2.0609271450, 1.8537379056, 4.3261222380, 0.3046098552, ... % V1 midget.LminusM 
38.3305740356, 0.5540646791, 0.0223900557, 0.1102570057, 0.8699638128, 0.6490083694, 0.7954258680, 0.6786286358, 2.9546277900, 3.8712201798, 1.6694139034, 0.4269341411, 0.4175314539, ... % V1 bistratified.S 
8.6200404167, 0.1425866365, 0.0761955976, 0.0038513422, 0.0244436264, 0.0003440857, 0.0001622677, 0.7500595918, 0.5690010263, 1.7807003553, 2.1076851417, 1.2945480328, 1.9776993246, ... % V1 parasol.LMS 
15.2464509010, 0.5579369307, 0.2844671965, 0.5465971231, 0.0005445719, 0.0013448000, 0.4263172865, 4.4351087595, 5.0054845618, 2.4345221910, 1.0904525195, 0.7295409434, 0.3000435808 ... % V1 midget.LMS 
 ]; 
        case 4
        otherwise
            % Fit the model with only RGC gain
            pMRI = fitMRIResponse(pMRI,stimulusDirections,studiedEccentricites,...
                studiedFreqs,v1Y,v1W,lgnY,lgnW,false,'v1');
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

%         % Show the cortical filter
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
%             a=gca; a.XTick = studiedFreqs; 
%             a.XTickLabel = arrayfun(@num2str, studiedFreqs, 'UniformOutput', 0);
%             a.XTickLabelRotation = 0;

        % Save the plot
        plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_MRIModelFit.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

    end

    figure
    for whichCell = 1:length(postReceptoralPaths)
        % Plot the surround suppression index vs. eccentricity
        subplot(1,2,1)
        paramIndices = 1+nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams: ...
            nUniqueParams+(whichCell-1)*(nFixedParams+nEccs*2)+nFixedParams+nEccs;
        plot(log10(studiedEccentricites),pMRI(paramIndices),['o' plotColor{whichCell}]);
        hold on
        for ee=1:nEccs
            plot([log10(studiedEccentricites(ee)) log10(studiedEccentricites(ee))],...
                [pMRI(paramIndices(ee))+pMRISEM(paramIndices(ee)) pMRI(paramIndices(ee))-pMRISEM(paramIndices(ee))],['-' plotColor{whichCell}]);
        end        
        plot(log10(studiedEccentricites),pMRI(paramIndices),['-' plotColor{whichCell}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Plot the gain index vs. eccentricity
        subplot(1,2,2)
        semilogy(log10(studiedEccentricites),pMRI(paramIndices+nEccs),['o' plotColor{whichCell}]);
        hold on
        for ee=1:nEccs
            plot([log10(studiedEccentricites(ee)) log10(studiedEccentricites(ee))],...
                [pMRI(paramIndices(ee)+nEccs)+pMRISEM(paramIndices(ee)+nEccs) pMRI(paramIndices(ee)+nEccs)-pMRISEM(paramIndices(ee)+nEccs)],['-' plotColor{whichCell}]);
        end        
        semilogy(log10(studiedEccentricites),pMRI(paramIndices+nEccs),['-' plotColor{whichCell}]);

        % Add the LGN gain values
        if whichCell < 4
            semilogy(0,pMRI(2+whichCell),['*' plotColor{whichCell}]);
        end

        xlim([-0.1 2]);
        xlabel('Eccentricity [log deg]');
        ylabel('Gain parameter');
    end

    % Save the plot
    plotName = [subjects{whichSub} '_MRIModelParams.pdf' ];
    saveas(gcf,fullfile(savePath,plotName));

end
