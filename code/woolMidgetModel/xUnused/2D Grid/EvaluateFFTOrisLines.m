%% EVALUATE FFT RESPONSES GIVEN SF LINES OR ORI CIRCLES
for t=1:length(thetas)
    OrientationLabel=strcat('Theta',num2str(thetas(t)),'deg');
    %Read the appropriate values from the RF FFTs
    lAmp=Lcpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    mAmp=Mcpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    normAmp=(max(max(Lcpower))+max(max(Mcpower)))/100;

    lPhase=Lcphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    mPhase=Mcphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);

    LcData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=lAmp/normAmp;
    McData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=mAmp/normAmp;

    LcData.(matlab.lang.makeValidName(OrientationLabel)).Phase=lPhase;%unwrap(lPhase);
    McData.(matlab.lang.makeValidName(OrientationLabel)).Phase=mPhase;%unwrap(mPhase);
end

for s=1:length(SpatFreqs)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    lAmp=Lcpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mAmp=Mcpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    normAmp=(max(max(Lcpower))+max(max(Mcpower)))/100;

    lPhase=Lcphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mPhase=Mcphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);

    LcData.(matlab.lang.makeValidName(SFLabel)).Amplitude=lAmp/normAmp;
    McData.(matlab.lang.makeValidName(SFLabel)).Amplitude=mAmp/normAmp;

    LcData.(matlab.lang.makeValidName(SFLabel)).Phase=lPhase;%unwrap(lPhase);
    McData.(matlab.lang.makeValidName(SFLabel)).Phase=mPhase;%unwrap(mPhase);
end

for t=1:length(thetas)
    OrientationLabel=strcat('Theta',num2str(thetas(t)),'deg');
    %Read the appropriate values from the RF FFTs
    lAmp=Lcpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    mAmp=Mcpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    normAmp=(max(max(Lspower))+max(max(Mspower)))/100;

    lPhase=Lsphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);
    mPhase=Msphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices);

    LsData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=lAmp/normAmp;
    MsData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=mAmp/normAmp;

    LsData.(matlab.lang.makeValidName(OrientationLabel)).Phase=lPhase;%unwrap(lPhase);
    MsData.(matlab.lang.makeValidName(OrientationLabel)).Phase=mPhase;%unwrap(mPhase);
end

for s=1:length(SpatFreqs)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    lAmp=Lspower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mAmp=Mspower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    normAmp=(max(max(Lspower))+max(max(Mspower)))/100;

    lPhase=Lsphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mPhase=Msphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);

    LsData.(matlab.lang.makeValidName(SFLabel)).Amplitude=lAmp/normAmp;
    MsData.(matlab.lang.makeValidName(SFLabel)).Amplitude=mAmp/normAmp;

    LsData.(matlab.lang.makeValidName(SFLabel)).Phase=lPhase;%unwrap(lPhase);
    MsData.(matlab.lang.makeValidName(SFLabel)).Phase=mPhase;%unwrap(mPhase);
end

%% PLOT 

set(0,'DefaultFigureWindowStyle','docked');
FiguresOpen=length(findobj(0,'type','figure'));

for f=1:length(thetas)
    OrientationLabel=strcat('Theta',num2str(thetas(f)),'deg');
    Freqs=Cell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).zFreqs;
    Line=Cell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).zIndices;
    PlotLc=LcData.(matlab.lang.makeValidName(OrientationLabel));
    PlotMc=McData.(matlab.lang.makeValidName(OrientationLabel));
    PlotLs=LsData.(matlab.lang.makeValidName(OrientationLabel));
    PlotMs=MsData.(matlab.lang.makeValidName(OrientationLabel));
    
    figure(f+1+FiguresOpen);
    centerstrength=strcat('Lc=',num2str(cpL*100),', Mc=',num2str(cpM*100));
    surroundstrength=strcat('Ls=',num2str(spL*100),', Ms=',num2str(spM*100));
    ori=strcat('Theta',num2str(thetas(f)),'deg');

    subplot(2,2,1);
    Ymax=max([PlotLc.Amplitude PlotMc.Amplitude]);
    Ymin=min([PlotLc.Amplitude PlotMc.Amplitude]);
    if Ymax==0
        Ymax=10;
    end
    PlotLAmp=loglog(Freqs,PlotLc.Amplitude,'ro');
    hold on;
    PlotMAmp=loglog(Freqs,PlotMc.Amplitude,'go');
    %PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Response amplitude (a.u.)'); 
    axis([0.005 15 .1 Ymax*10]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    legend([PlotLAmp,PlotMAmp],'Location','southwest')
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
    %     set(PlotSAmp,...
    %         'DisplayName','S',...
    %         'LineWidth',.5,...
    %         'LineStyle',':',...
    %         'Color','k',...
    %         'MarkerFaceColor',[0 0.45 0.74])

    subplot(2,2,3);
    PlotLPhase=semilogx(Freqs,abs(PlotLc.Phase),'ro');
    hold on;
    PlotMPhase=semilogx(Freqs,abs(PlotMc.Phase),'go');
    %PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Phase (degrees)'); 
    axis([0.005 15 -20 200]);
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
    %     set(PlotSPhase,...
    %         'DisplayName','S',...
    %         'LineWidth',.5,...
    %         'LineStyle',':',...
    %         'Color','k',...
    %         'MarkerFaceColor',[0 0.45 0.74])

    subplot(2,2,2);
    Ymax=max([PlotLs.Amplitude PlotMs.Amplitude]);
    Ymin=min([PlotLs.Amplitude PlotMs.Amplitude]);
    if Ymax==0
        Ymax=10;
    end
    PlotLAmp=loglog(Freqs,PlotLs.Amplitude,'ro');
    hold on;
    PlotMAmp=loglog(Freqs,PlotMs.Amplitude,'go');
    %PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Response amplitude (a.u.)'); 
    axis([0.005 15 .1 Ymax*10]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    legend([PlotLAmp,PlotMAmp],'Location','southwest')
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
    %     set(PlotSAmp,...
    %         'DisplayName','S',...
    %         'LineWidth',.5,...
    %         'LineStyle',':',...
    %         'Color','k',...
    %         'MarkerFaceColor',[0 0.45 0.74])

    subplot(2,2,4);
    PlotLPhase=semilogx(Freqs,abs(PlotLs.Phase),'ro');
    hold on;
    PlotMPhase=semilogx(Freqs,abs(PlotMs.Phase),'go');
    %PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Phase (degrees)'); 
    axis([0.005 15 -20 200]);
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
    %     set(PlotSPhase,...
    %         'DisplayName','S',...
    %         'LineWidth',.5,...
    %         'LineStyle',':',...
    %         'Color','k',...
    %         'MarkerFaceColor',[0 0.45 0.74])

    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,{centerstrength;surroundstrength;ori},'HorizontalAlignment','center','VerticalAlignment', 'top');
end
set(0,'DefaultFigureWindowStyle','default');