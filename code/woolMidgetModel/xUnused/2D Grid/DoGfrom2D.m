function Cell=DoGfrom2D(Em)
%% IMPORT DATA
Cell=GenerateCell(Em);

cpAll=Cell.CenterWeights; %proportion of L and M cones
spAll=Cell.SurroundWeights;

%% MODEL UNIVERSE

Fs=512; %Sampling rate
numdeg=512; %degrees
x=-numdeg/2:1/Fs:numdeg/2-(1/Fs); %length(x) = numdeg * Fs
%Em=5; %Eccentricity in mm
Ed=(Em*1000)/200; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)
w=length(x);
Freq=(-w/2:(w/2)-1)*(Fs/w);

SigmasToRFRadius=1; %Let's say that ~3 SDs make up a RF radius
CSRadRatio=6; %Let's also say that a surround radius is ~6x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Finally, let's say that surround gain is ~3/4 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.3736*Em^2 + 0.01126*Em + 1.026);
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

DFRadius=0.002738*Em.^1.327; %DF radius data, in mm
DFRadiusDeg=((0.002738*Em.^1.327)*1000)/200; %This is in degrees.
CenterSigma=DFRadiusDeg/SigmasToRFRadius; %Here we assign a sigma for our cell given a particular eccentricity.
SurroundSigma=CSRadRatio*CenterSigma; %Ditto.

%ON or OFF cell? -1: OFF, +1: ON
CellType=1;

%% RF STRUCTURE

%RF center (c) and surround (s) gains
crL=CellType*100*cpAll(1);
crM=CellType*100*cpAll(2);
crS=CellType*100*cpAll(3);

srL=-CellType*100*spAll(1)*CSStrengthRatio;
srM=-CellType*100*spAll(2)*CSStrengthRatio;
srS=-CellType*100*spAll(3)*CSStrengthRatio;

% RF center (c) and surround (s) sizes
LcSigma=CenterSigma; %Center radius (degrees)
LcAmp=crL; %Strength of sensitivity
Lc=LcAmp*normpdf(x,0,LcSigma);

McSigma=CenterSigma;
McAmp=crM;
Mc=McAmp*normpdf(x,0,McSigma);

ScSigma=0;
ScAmp=crS;
Sc=ScAmp*normpdf(x,0,ScSigma);

LsSigma=SurroundSigma;
LsAmp=srL;
Ls=LsAmp*normpdf(x,0,LsSigma);

MsSigma=SurroundSigma;
MsAmp=srM;
Ms=MsAmp*normpdf(x,0,MsSigma);

SsSigma=0;
SsAmp=crS;
Ss=SsAmp*normpdf(x,0,SsSigma);

% Total cell response:
CellResp=Lc+Ls+Mc+Ms+Sc+Ss;

%% PRESENT THE STIMULUS

Freqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32 64 128];
%SFLabels={'SF=1/128' 'SF=1/64' 'SF=1/32' 'SF=1/16' 'SF=1/8' 'SF=1/4' 'SF=1/2' 'SF=1' 'SF=2' 'SF=4' 'SF=8' 'SF=16' 'SF=32' 'SF=64' 'SF=128'};

for t=1:length(Freqs)
    SpatFreq=Freqs(t); %Spatial frequency (cpd) **THIS MUST BE A POWER OF 2**
    T=1/SpatFreq; %Fundamental period
    TotalPeriods=(max(x)-min(x))/T;
    phi=0;
    A=1;
    L_Iso=A*cos(2*pi*SpatFreq*x+phi);
    M_Iso=A*cos(2*pi*SpatFreq*x+phi);
    %S_Iso=A*cos(2*pi*SpatFreq*x+phi);

    %%%%%%%%%% OUTPUT RESPONSE %%%%%%%%%%%%%

%     LcResp=1/Fs*conv(Lc,L_Iso,'same');
%     McResp=1/Fs*conv(Mc,M_Iso,'same');
%     %ScResp=1/Fs*conv(Sc,S_Iso,'same');
% 
%     LsResp=1/Fs*conv(Ls,L_Iso,'same');
%     MsResp=1/Fs*conv(Ms,M_Iso,'same');
%     %SsResp=1/Fs*conv(Ss,S_Iso,'same');
    
    LcResp=1/Fs*ifft(fft(Lc).*fft(L_Iso));
    McResp=1/Fs*ifft(fft(Mc).*fft(M_Iso));
    LsResp=1/Fs*ifft(fft(Ls).*fft(L_Iso));
    MsResp=1/Fs*ifft(fft(Ms).*fft(M_Iso));

    TotalLResp=LcResp+LsResp;
    TotalMResp=McResp+MsResp;
    %TotalSResp=ScResp+SsResp;
    %TotalResp=TotalLResp+TotalMResp+TotalSResp;
    TotalResp=TotalLResp+TotalMResp;
    LMDiff=TotalLResp-TotalMResp;

    bound=ceil(max([LsSigma MsSigma SsSigma])*3*SpatFreq)/SpatFreq; %Determine the gaussian width, round up to nearest whole cycle
    window=-(numdeg/2)+bound:1/Fs:(numdeg/2)-bound-1/Fs;
    frame=(length(x)-length(window))/2;
    wx=x(frame+1:(end-frame));
    ww=length(wx);
    wFreq=(-ww/2:(ww/2)-1)*(Fs/w);
    
%     frame=0;
%     ww=length(x);
%     wFreq=Freq;

    wLcResp=LcResp(frame+1:(end-frame));
    wMcResp=McResp(frame+1:(end-frame));
    %wScResp=ScResp(frame+1:(end-frame));

    wLsResp=LsResp(frame+1:(end-frame));
    wMsResp=MsResp(frame+1:(end-frame));
    %wSsResp=SsResp(frame+1:(end-frame));

    wTotalLResp=TotalLResp(frame+1:(end-frame));
    wTotalMResp=TotalMResp(frame+1:(end-frame));
    %wTotalSResp=TotalSResp(frame+1:(end-frame));
    wTotalResp=TotalResp(frame+1:(end-frame));

    wLMDiff=LMDiff(frame+1:(end-frame));

    %%%%%%%%%% FFT ANALYSIS %%%%%%%%%%%%%%%%

    TotalXf=(fft(wTotalResp,ww));
    TotalPowerXf=abs(TotalXf)/floor(ww/2);
    TotalPhaseXf= angle(TotalXf)*180/pi;
    [Amp, Index] = max(TotalPowerXf);
    TotalData.Amplitude(t) = Amp;
    TotalData.Frequency(t) = numdeg/2 + Freq(Index);
    TotalData.Phase(t) = TotalPhaseXf(Index);

    LMDiffXf=(fft(wLMDiff,ww));
    LMDiffPowerXf=abs(LMDiffXf)/floor(ww/2);
    LMDiffPhaseXf= angle(LMDiffXf)*180/pi;
    [Amp, Index] = max(LMDiffPowerXf);
    LMDiffData.Amplitude(t) = Amp;
    LMDiffData.Frequency(t) = numdeg/2 + Freq(Index);
    LMDiffData.Phase(t) = LMDiffPhaseXf(Index);

    TotalLXf=(fft(wTotalLResp,ww));
    TotalLPowerXf=abs(TotalLXf)/floor(ww/2);
    TotalLPhaseXf= angle(TotalLXf)*180/pi;
    [Amp, Index] = max(TotalLPowerXf);
    TotalLData.Amplitude(t) = Amp;
    TotalLData.Frequency(t) = numdeg/2 + Freq(Index);
    TotalLData.Phase(t) = TotalLPhaseXf(Index);

    TotalMXf=(fft(wTotalMResp,ww));
    TotalMPowerXf=abs(TotalMXf)/floor(ww/2);
    TotalMPhaseXf= angle(TotalMXf)*180/pi;
    [Amp, Index] = max(TotalMPowerXf);
    TotalMData.Amplitude(t) = Amp;
    TotalMData.Frequency(t) = numdeg/2 + Freq(Index);
    TotalMData.Phase(t) = TotalMPhaseXf(Index);

%         TotalSXf=(fft(wTotalSResp,ww));
%         TotalSPowerXf=abs(TotalSXf)/floor(ww/2);
%         TotalSPhaseXf= angle(TotalSXf)*180/pi;
%         [Amp, Index] = max(TotalSPowerXf);
%         TotalSData.Amplitude(t) = Amp;
%         TotalSData.Frequency(t) = numdeg/2 + Freq(Index);
%         TotalSData.Phase(t) = TotalSPhaseXf(Index);

    LcXf=(fft(LcResp,w));
    LcPowerXf=abs(LcXf)/floor(w/2);
    LcPhaseXf= angle(LcXf)*180/pi;
    [Amp, Index] = max(LcPowerXf);
    LcData.Amplitude(t) = Amp;
    LcData.Frequency(t) = numdeg/2 + Freq(Index);
    LcData.Phase(t) = LcPhaseXf(Index);

    McXf=(fft(McResp,w));
    McPowerXf=abs(McXf)/floor(w/2);
    McPhaseXf= angle(McXf)*180/pi;
    [Amp, Index] = max(McPowerXf);
    McData.Amplitude(t) = Amp;
    McData.Frequency(t) = numdeg/2 + Freq(Index);
    McData.Phase(t) = McPhaseXf(Index);

%         ScXf=(fft(ScResp,w));
%         ScPowerXf=abs(ScXf)/floor(w/2);
%         ScPhaseXf= angle(ScXf)*180/pi;
%         [Amp, Index] = max(ScPowerXf);
%         ScData.Amplitude(t) = Amp;
%         ScData.Frequency(t) = numdeg/2 + Freq(Index);
%         ScData.Phase(t) = ScPhaseXf(Index);

    LsXf=(fft(LsResp,w));
    LsPowerXf=abs(LsXf)/floor(w/2);
    LsPhaseXf= angle(LsXf)*180/pi;
    [Amp, Index] = max(LsPowerXf);
    LsData.Amplitude(t) = Amp;
    LsData.Frequency(t) = numdeg/2 + Freq(Index);
    LsData.Phase(t) = LsPhaseXf(Index);

    MsXf=(fft(MsResp,w));
    MsPowerXf=abs(MsXf)/floor(w/2);
    MsPhaseXf= angle(MsXf)*180/pi;
    [Amp, Index] = max(MsPowerXf);
    MsData.Amplitude(t) = Amp;
    MsData.Frequency(t) = numdeg/2 + Freq(Index);
    MsData.Phase(t) = MsPhaseXf(Index);

%         SsXf=(fft(SsResp,w));
%         SsPowerXf=abs(SsXf)/floor(w/2);
%         SsPhaseXf= angle(SsXf)*180/pi;
%         [Amp, Index] = max(SsPowerXf);
%         SsData.Amplitude(t) = Amp;
%         SsData.Frequency(t) = numdeg/2 + Freq(Index);
%         SsData.Phase(t) = SsPhaseXf(Index);

% %Plot the output response (total and component)
% figure;
% hold on;
% TotalPlot=plot(wx,wTotalResp,'DisplayName','TotalPlot');
% LcPlot=plot(wx,wLcResp,'DisplayName','LcPlot');
% LsPlot=plot(wx,wLsResp,'DisplayName','LsPlot');
% McPlot=plot(wx,wMcResp,'DisplayName','McPlot');
% MsPlot=plot(wx,wMsResp,'DisplayName','MsPlot');
% % ScPlot=plot(wx,wScResp,'DisplayName','ScPlot');
% % SsPlot=plot(wx,wSsResp,'DisplayName','SsPlot');

end



%% TEST FOR CHROMATIC OPPONENCY 

%Assess chromatic opponency using the phase relationship
AmpIndex=TotalData.Amplitude(1)-LMDiffData.Amplitude(1);
PhaseIndex=abs(TotalLData.Phase(1))-abs(TotalMData.Phase(1));

if AmpIndex<0 %L-M response is greater than L+M response at full-field stimulation
    if PhaseIndex<-175 %L phase = 0, M phase = 180
        ChromTag='L-dominated';
    elseif PhaseIndex>175 %M phase = 0, L phase = 180
        ChromTag='M-dominated';
    end
elseif AmpIndex>0 %L+M response is greater than/equal to L-M response at full-field stimulation
    ChromTag='Achromatic';
end

%% SAVE DATA

Cell.ChromTag=ChromTag;

Cell.Convolution.DFRadiusDeg=DFRadiusDeg;
Cell.Convolution.SpatialFrequencies=Freqs;

Cell.Convolution.LResponse.CenterSigma=LcSigma;
Cell.Convolution.LResponse.CenterAmplitude=crL;
Cell.Convolution.LResponse.SurroundSigma=LsSigma;
Cell.Convolution.LResponse.SurroundAmplitude=srL;

Cell.Convolution.MResponse.CenterSigma=McSigma;
Cell.Convolution.MResponse.CenterAmplitude=crM;
Cell.Convolution.MResponse.SurroundSigma=MsSigma;
Cell.Convolution.MResponse.SurroundAmplitude=srM;

Cell.Convolution.SResponse.CenterSigma=ScSigma;
Cell.Convolution.SResponse.CenterAmplitude=crS;
Cell.Convolution.SResponse.SurroundSigma=SsSigma;
Cell.Convolution.SResponse.SurroundAmplitude=srS;

Cell.ResponseFunctions.LIsoResponse.Amplitude=TotalLData.Amplitude;
Cell.ResponseFunctions.LIsoResponse.Phase=TotalLData.Phase;
Cell.ResponseFunctions.MIsoResponse.Amplitude=TotalMData.Amplitude;
Cell.ResponseFunctions.MIsoResponse.Phase=TotalMData.Phase;
Cell.ResponseFunctions.LMSumResponse.Amplitude=TotalData.Amplitude;
Cell.ResponseFunctions.LMSumResponse.Phase=TotalData.Phase;
Cell.ResponseFunctions.LMDiffResponse.Amplitude=LMDiffData.Amplitude;
Cell.ResponseFunctions.LMDiffResponse.Phase=LMDiffData.Phase;

%% PLOTTING

% %RF FIELD COMPONENTS
% figure;
% hold on; 
% TotPlot=plot(x,CellResp,'DisplayName','Total response');
% LcPlot=plot(x,Lc,'DisplayName','L-cone center response'); 
% LsPlot=plot(x,Ls,'DisplayName','L-cone surround response'); 
% McPlot=plot(x,Mc,'DisplayName','M-cone center response'); 
% MsPlot=plot(x,Ms,'DisplayName','M-cone surround response');
% %ScPlot=plot(x,Sc,'DisplayName','S-cone center response'); 
% %SsPlot=plot(x,Ss,'DisplayName','S-cone surround response');
% axis([-2 2 min([Ls Ms])*1.5 max([Lc Mc])*1.05]);
% legend([LcPlot,McPlot,LsPlot,MsPlot]);

% %AMPLITUDE AND PHASE OUTPUTS FOR ISO, COMBO, DIFF
% figure;
% centerstrength=strcat('Lc=',num2str(crL),', Mc=',num2str(crM),', Sc=',num2str(crS));
% surroundstrength=strcat('Ls=',num2str(srL),', Ms=',num2str(srM),', Ss=',num2str(srS));
% 
% subplot(2,3,1);
% Ymax=max([TotalLData.Amplitude TotalMData.Amplitude]);
% if Ymax==0
%     Ymax=10;
% end
% PlotLAmp=loglog(Freqs,TotalLData.Amplitude,'ro');
% hold on;
% PlotMAmp=loglog(Freqs,TotalMData.Amplitude,'go');
% %PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Response amplitude (a.u.)'); 
% axis([0.005 150 1 100]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% legend([PlotLAmp,PlotMAmp],'Location','southwest')
% set(PlotLAmp,...
%     'DisplayName','L',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[204 0 0]/255)
% set(PlotMAmp,...
%     'DisplayName','M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[119 172 48]/255)
% %     set(PlotSAmp,...
% %         'DisplayName','S',...
% %         'LineWidth',.5,...
% %         'LineStyle',':',...
% %         'Color','k',...
% %         'MarkerFaceColor',[0 0.45 0.74])
% 
% subplot(2,3,4);
% PlotLPhase=semilogx(Freqs,abs(TotalLData.Phase),'ro');
% hold on;
% PlotMPhase=semilogx(Freqs,abs(TotalMData.Phase),'go');
% %PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Phase (degrees)'); 
% axis([0.005 150 -20 200]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% %legend([PlotLPhase,PlotMPhase],'Location','northeast')
% set(PlotLPhase,...
%     'DisplayName','L',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[204 0 0]/255)
% set(PlotMPhase,...
%     'DisplayName','M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[119 172 48]/255)
% %     set(PlotSPhase,...
% %         'DisplayName','S',...
% %         'LineWidth',.5,...
% %         'LineStyle',':',...
% %         'Color','k',...
% %         'MarkerFaceColor',[0 0.45 0.74])
% 
% subplot(2,3,2);
% Ymax=max([TotalData.Amplitude]);
% if Ymax==0
%     Ymax=10;
% end
% PlotTotAmp=loglog(Freqs,TotalData.Amplitude,'ro');
% hold on;
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Response amplitude (a.u.)'); 
% axis([0.005 150 1 100]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% legend(PlotTotAmp,'Location','southwest')
% set(PlotTotAmp,...
%     'DisplayName','L+M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[.8 .8 .8])
% 
% subplot(2,3,5);
% PlotTotPhase=semilogx(Freqs,abs(TotalData.Phase),'ro');
% hold on;
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Phase (degrees)'); 
% axis([0.005 150 -20 200]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% %legend(PlotTotPhase,'Location','northeast')
% set(PlotTotPhase,...
%     'DisplayName','L+M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[.8 .8 .8])
% 
% subplot(2,3,3);
% Ymax=max([LMDiffData.Amplitude]);
% if Ymax==0
%     Ymax=10;
% end
% PlotLMDiffAmp=loglog(Freqs,LMDiffData.Amplitude,'ro');
% hold on;
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Response amplitude (a.u.)'); 
% axis([0.005 150 1 100]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% legend(PlotLMDiffAmp,'Location','southwest')
% set(PlotLMDiffAmp,...
%     'DisplayName','L-M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[.31 .31 .31])
% 
% subplot(2,3,6);
% PlotLMDiffPhase=semilogx(Freqs,abs(LMDiffData.Phase),'ro');
% hold on;
% %title({centerstrength;surroundstrength});
% xlabel('Stimulus frequency (cpd)');
% ylabel('Phase (degrees)'); 
% axis([0.005 150 -20 200]);
% set(gca,'TickDir','in','TickLength', [.005 .005]);box off
% %legend(PlotTotPhase,'Location','northeast')
% set(PlotLMDiffPhase,...
%     'DisplayName','L-M',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[.31 .31 .31])
% 
% axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 1,{centerstrength;surroundstrength},'HorizontalAlignment','center','VerticalAlignment', 'top');


%% WORKSHOP

%Dendritic tree size (degrees)

% Em=1; %Eccentricity in mm
% Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)
% 
% %Center cone count
% ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
% ConesToSurround=6*ConesToCenter;
% 
% pL=0.6;
% pM=0.3;
% pS=0.1;
% 
% cnL=CellType*sum(binornd(1:1:1,pL,1,ConesToCenter));
% cnM=CellType*sum(binornd(1:1:1,pM,1,ConesToCenter));
% cnS=CellType*sum(binornd(1:1:1,pS,1,ConesToCenter));
% 
% snL=-CellType*sum(binornd(1:1:1,pL,1,ConesToSurround));
% snM=-CellType*sum(binornd(1:1:1,pM,1,ConesToSurround));
% snS=-CellType*sum(binornd(1:1:1,pS,1,ConesToSurround));
% 
% disp(strcat(['Cones to center: ',num2str(ConesToCenter),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
% disp(strcat(['Cones to surround: ',num2str(ConesToSurround),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));

% %Plot the output response (total and component)
% figure(2);
% hold on;
% TotalPlot=plot(wx,wTotalResp,'DisplayName','TotalPlot');
% LcPlot=plot(wx,wLcResp,'DisplayName','LcPlot');
% LsPlot=plot(wx,wLsResp,'DisplayName','LsPlot');
% McPlot=plot(wx,wMcResp,'DisplayName','McPlot');
% MsPlot=plot(wx,wMsResp,'DisplayName','MsPlot');
% ScPlot=plot(wx,wScResp,'DisplayName','ScPlot');
% SsPlot=plot(wx,wSsResp,'DisplayName','SsPlot');
