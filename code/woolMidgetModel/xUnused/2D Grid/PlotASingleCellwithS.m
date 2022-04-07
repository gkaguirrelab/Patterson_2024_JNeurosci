PlotCell=Cell;
CellLAmpData=PlotCell.ResponseFunctions.LIsoResponse.Theta0deg.Amplitude;
CellMAmpData=PlotCell.ResponseFunctions.MIsoResponse.Theta0deg.Amplitude;
CellSAmpData=PlotCell.ResponseFunctions.SIsoResponse.Theta0deg.Amplitude;
CellLPhaseData=PlotCell.ResponseFunctions.LIsoResponse.Theta0deg.Phase;
CellMPhaseData=PlotCell.ResponseFunctions.MIsoResponse.Theta0deg.Phase;
CellSPhaseData=PlotCell.ResponseFunctions.SIsoResponse.Theta0deg.Phase;
CellLMSumAmpData=PlotCell.ResponseFunctions.LMSumResponse.Theta0deg.Amplitude;
CellLMSumPhaseData=PlotCell.ResponseFunctions.LMSumResponse.Theta0deg.Phase;
CellLMDiffAmpData=PlotCell.ResponseFunctions.LMDiffResponse.Theta0deg.Amplitude;
CellLMDiffPhaseData=PlotCell.ResponseFunctions.LMDiffResponse.Theta0deg.Phase;

CellCenterWeight=PlotCell.CenterWeights;
CellSurroundWeight=PlotCell.SurroundWeights;

Freqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32 64 128];
SFLabels={'SF=1/128' 'SF=1/64' 'SF=1/32' 'SF=1/16' 'SF=1/8' 'SF=1/4' 'SF=1/2' 'SF=1' 'SF=2' 'SF=4' 'SF=8' 'SF=16' 'SF=32' 'SF=64' 'SF=128'};
centerstrength=strcat('Lc=',num2str(CellCenterWeight(1)),', Mc=',num2str(CellCenterWeight(2)),', Sc=',num2str(CellCenterWeight(3)));
surroundstrength=strcat('Ls=',num2str(CellSurroundWeight(1)),', Ms=',num2str(CellSurroundWeight(2)),', Ss=',num2str(CellSurroundWeight(3)));

figure;
subplot(2,3,1);
Ymax=max([CellLAmpData CellMAmpData]);
if Ymax==0
    Ymax=10;
end
PlotLAmp=loglog(Freqs,CellLAmpData,'ro');
hold on;
PlotMAmp=loglog(Freqs,CellMAmpData,'go');
PlotSAmp=loglog(Freqs,CellSAmpData,'bo');
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend([PlotLAmp,PlotMAmp,PlotSAmp],'Location','southwest')
set(PlotLAmp,...
    'DisplayName','L',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[204 0 0]/255)
set(PlotMAmp,...
    'DisplayName','M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[119 172 48]/255)
set(PlotSAmp,...
    'DisplayName','S',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[0 0.45 0.74])
if CellSAmpData(1)==0
    set(PlotSAmp,'visible', 'off')
end

subplot(2,3,4);
PlotLPhase=semilogx(Freqs,abs(CellLPhaseData),'ro');
hold on;
PlotMPhase=semilogx(Freqs,abs(CellMPhaseData),'go');
PlotSPhase=semilogx(Freqs,abs(CellSPhaseData),'bo');
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend([PlotLPhase,PlotMPhase],'Location','northeast')
set(PlotLPhase,...
    'DisplayName','L',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[204 0 0]/255)
set(PlotMPhase,...
    'DisplayName','M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[119 172 48]/255)
set(PlotSPhase,...
    'DisplayName','S',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[0 0.45 0.74])
if CellSAmpData(1)==0
    set(PlotSPhase,'visible', 'off')
end

subplot(2,3,2);
Ymax=max([CellLMSumAmpData]);
if Ymax==0
    Ymax=10;
end
PlotTotAmp=loglog(Freqs,CellLMSumAmpData,'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend(PlotTotAmp,'Location','southwest')
set(PlotTotAmp,...
    'DisplayName','L+M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.8 .8 .8])

subplot(2,3,5);
PlotTotPhase=semilogx(Freqs,abs(CellLMSumPhaseData),'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend(PlotTotPhase,'Location','northeast')
set(PlotTotPhase,...
    'DisplayName','L+M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.8 .8 .8])

subplot(2,3,3);
Ymax=max([CellLMDiffAmpData]);
if Ymax==0
    Ymax=10;
end
PlotLMDiffAmp=loglog(Freqs,CellLMDiffAmpData,'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 1 100]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
legend(PlotLMDiffAmp,'Location','southwest')

set(PlotLMDiffAmp,...
    'DisplayName','L-M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.31 .31 .31])

subplot(2,3,6);
PlotLMDiffPhase=semilogx(Freqs,abs(CellLMDiffPhaseData),'ro');
hold on;
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Phase (degrees)'); 
axis([0.005 150 -20 200]);
set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%legend(PlotTotPhase,'Location','northeast')
set(PlotLMDiffPhase,...
    'DisplayName','L-M',...
    'LineWidth',.5,...
    'LineStyle',':',...
    'Color','k',...
    'MarkerFaceColor',[.31 .31 .31])

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,{centerstrength;surroundstrength},'HorizontalAlignment','center','VerticalAlignment', 'top');