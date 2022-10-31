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

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots');

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
plotColor = mriFullResultSet.meta.plotColor;
nFixedParams = mriFullResultSet.meta.nFixedParams;
nFloatByEccParams = mriFullResultSet.meta.nFloatByEccParams;
nUniqueParams = mriFullResultSet.meta.nUniqueParams;
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
freqsForPlotting = logspace(0,2,50);
nFreqsForPlotting = length(freqsForPlotting);
cellClassOrder = {'midgetChrom','parasol','bistratified','midgetAchrom'};
chromAchromIndex = [1 1 2];
nParamsPerCellBlock = nFixedParams+nEccs*2;

for whichSub = 1:length(subjects)

    %load(fullfile(loadPath,subjects{whichSub},['mriTemporalModel_' regexprep(num2str(1:12),' +','-') '.mat']),'mriTemporalModel');

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI);
    v1Y = mriFullResultSet.(subjects{whichSub}).v1YMean;
    v1Ylow = mriFullResultSet.(subjects{whichSub}).v1Y_lowCI;
    v1Yhigh = mriFullResultSet.(subjects{whichSub}).v1Y_highCI;
    lgnY = mriFullResultSet.(subjects{whichSub}).lgnYMean;
    lgnYlow = mriFullResultSet.(subjects{whichSub}).lgnY_lowCI;
    lgnYhigh = mriFullResultSet.(subjects{whichSub}).lgnY_highCI;

    [~,v1YFitMatrix] = assembleV1ResponseAcrossStimsAndEcc(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel,nUniqueParams,nFixedParams);
    [~,lgnYFitMatrix] = assembleLGNResponseAcrossStims(pMRI,stimulusDirections,freqsForPlotting,rgcTemporalModel);

    for whichStim = 1:length(stimulusDirections)
        figure

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
        ylim([-0.5 4]);

        % Show the cortical filter
        subplot(2,4,4)
        secondOrderFc = pMRI(5+(chromAchromIndex(whichStim)-1)*3+1);
        secondOrderQ = pMRI(5+(chromAchromIndex(whichStim)-1)*3+2);

        syms f; rf = stageSecondOrderLP(f,secondOrderFc,secondOrderQ);
        myFreqs = logspace(log10(0.5),log10(100),101);
        ttfComplex = double(subs(rf,myFreqs));
        gainVals = abs(ttfComplex);
        semilogx(myFreqs,gainVals,['-.' plotColor{whichStim}]);
        xlabel('frequency [Hz]'); ylabel('gain');

        % Save the plot
        plotName = [stimulusDirections{whichStim} '_' subjects{whichSub} '_MRIModelFit.pdf' ];
        saveas(gcf,fullfile(savePath,plotName));

    end

    figure
    for whichCell = 1:length(cellClassOrder)
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
