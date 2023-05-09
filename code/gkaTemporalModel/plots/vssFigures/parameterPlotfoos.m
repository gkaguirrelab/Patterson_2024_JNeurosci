% scriptCreatePlots

% Housekeeping
clear
close all

% Where will we save the plots
savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots',modelType);



subjectLineSpec = {'-','--'};
subjectSymbol = {'.','o'};

paramValsFig = figure();
figuresize(800,600,'pt');

% Loop over subjects and make plots
for whichSub = 1:length(subjects)

    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

    % Make the paramVals figure active
    figure(paramValsFig)

    % Plot the surround suppression index
    subOrder = [2,3,1];
    for whichStim = 1:length(stimulusDirections)

        subplot(2,4,subOrder(whichStim));

        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed;
        v1SurroundIndex = pMRI(startIdx+1:startIdx+nEccs);
        v1SurroundIndexSEM = pMRISEM(startIdx+1:startIdx+nEccs);

        % Plot the surround suppression index vs. eccentricity
        plot(log10(studiedEccentricites),v1SurroundIndex,[subjectLineSpec{whichSub} plotColor{whichStim}],'LineWidth',2);
        hold on

        % Add error bars
        X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
        Y = [v1SurroundIndex-v1SurroundIndexSEM, fliplr(v1SurroundIndex+v1SurroundIndexSEM)];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.2);

        % Clean up
        semilogy([-0.5 2],[0 0],':k');
        a=gca;
        a.XTick = log10([2,4,8,16,32,54]);
        a.XTickLabel = {'2','4','8','16','32','64'};
        xlabel('Eccentricity [deg]');
        ylabel('Suppression index');
        xlim([0 2]);
        ylim([-1 1]);
        box off
    end

    % Calculate the gain ratios
    v1GainVals = {}; lgnGainVals = {};
    for whichStim = 1:length(stimulusDirections)
        startIdx = paramCounts.unique + paramCounts.lgn*nCells + (whichStim-1)*paramCounts.v1total + paramCounts.v1fixed + nEccs + 1;
        v1GainVals(whichStim) = {mriFullResultSet.(subjects{whichSub}).pMRI(:,startIdx:startIdx+nEccs-1)};
        startIdx = paramCounts.unique + whichStim*paramCounts.lgn;
        lgnGainVals(whichStim) = {mriFullResultSet.(subjects{whichSub}).pMRI(:,startIdx)};
    end

    % Plot the luminance gain ratio relative to fovea
    for whichStim = 1:length(stimulusDirections)
        subplot(2,4,4+subOrder(whichStim))
        meanV1Val = mean(v1GainVals{whichStim});
        semV1Val = std(v1GainVals{whichStim});
        semilogy(log10(studiedEccentricites),meanV1Val,[subjectLineSpec{whichSub} plotColor{whichStim}],'LineWidth',2);
        hold on
        X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
        Y = [meanV1Val+semV1Val, fliplr(meanV1Val-semV1Val)];
        p = patch(X,Y,plotColor{whichStim});
        set(p,'edgecolor','none','facealpha',0.2);
        semilogy([-0.5 2],[1 1],':k');
        xlim([0 2]);
        a=gca;
        a.XTick = log10([2,4,8,16,32,54]);
        a.XTickLabel = {'2','4','8','16','32','64'};
        ylim([10^-2 10^2]);
        xlabel('Eccentricity [deg]');
        ylabel('Retinal gain');
        box off
    end

    % Plot the relative gain ratios
    subplot(2,4,8)

    % achromatic vs. achromatic
    refStim = 3;
    meanV1Val = mean(v1GainVals{refStim}./v1GainVals{refStim});
    semV1Val = std(v1GainVals{refStim}./mean(v1GainVals{refStim}));
    meanLGNVal = mean(lgnGainVals{refStim}./lgnGainVals{refStim});
    semLGNVal = std(lgnGainVals{refStim}./mean(lgnGainVals{refStim}));
    semilogy(log10(studiedEccentricites),meanV1Val,[subjectLineSpec{whichSub} 'k'],'LineWidth',2);
    hold on

    % Omit the LGN point if the value does not differ from zero
    if (meanLGNVal-semLGNVal)>0 && showLGNValsOnGainPlot
        semilogy(0,meanLGNVal,[subjectSymbol{whichSub} 'k'])
        semilogy([0 0],[meanLGNVal-semLGNVal,meanLGNVal+semLGNVal,'-k'])
    end
    X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
    Y = [meanV1Val+semV1Val, fliplr(meanV1Val-semV1Val)];
    p = patch(X,Y,'k');
    set(p,'edgecolor','none','facealpha',0.2);

    % RG chrom vs. achromatic
    meanV1Val = mean(v1GainVals{1}./v1GainVals{refStim});
    semV1Val = std(v1GainVals{1}./v1GainVals{refStim});
    meanLGNVal = mean(lgnGainVals{1}./lgnGainVals{refStim});
    semLGNVal = std(lgnGainVals{1}./lgnGainVals{refStim});
    semilogy(log10(studiedEccentricites),meanV1Val,[subjectLineSpec{whichSub} 'r'],'LineWidth',2);

    % Omit the LGN point if the value does not differ from zero
    if (meanLGNVal-semLGNVal)>0 && showLGNValsOnGainPlot
        semilogy(-0.1,meanLGNVal,[subjectSymbol{whichSub} 'r'])
        semilogy([-0.1 -0.1],[meanLGNVal-semLGNVal,meanLGNVal+semLGNVal,'-r'])
    end
    X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
    Y = [meanV1Val+semV1Val, fliplr(meanV1Val-semV1Val)];
    p = patch(X,Y,'r');
    set(p,'edgecolor','none','facealpha',0.2);

    % BY chrom vs. achromatic
    meanV1Val = mean(v1GainVals{2}./v1GainVals{refStim});
    semV1Val = std(v1GainVals{2}./v1GainVals{refStim});
    meanLGNVal = mean(lgnGainVals{2}./lgnGainVals{refStim});
    semLGNVal = std(lgnGainVals{2}./lgnGainVals{refStim});
    semilogy(log10(studiedEccentricites),meanV1Val,[subjectLineSpec{whichSub} 'b'],'LineWidth',2);

    % Omit the LGN point if the value does not differ from zero
    if (meanLGNVal-semLGNVal)>0 && showLGNValsOnGainPlot
        semilogy(0.1,meanLGNVal,[subjectSymbol{whichSub} 'b'])
        semilogy([0.1 0.1],[meanLGNVal-semLGNVal,meanLGNVal+semLGNVal],'-b')
    end
    X = [log10(studiedEccentricites) fliplr(log10(studiedEccentricites))];
    Y = [meanV1Val+semV1Val, fliplr(meanV1Val-semV1Val)];
    p = patch(X,Y,'b');
    set(p,'edgecolor','none','facealpha',0.2);

    % Clean up
    semilogy([-0.5 2],[1 1],':k');
    xlim([0 2]);
    ylim([10^-1 10^1.5]);
    a=gca;
    a.XTick = log10([2,4,8,16,32,54]);
    a.XTickLabel = {'2','4','8','16','32','64'};
    xlabel('Eccentricity [deg]');
    ylabel('Retinal gain relative to luminance');
    box off


end


% Save across-subject parameter plot
plotName = [paramSearch '_RelativeGainAcrossSubjects.pdf' ];
saveas(paramValsFig,fullfile(savePath,plotName));

