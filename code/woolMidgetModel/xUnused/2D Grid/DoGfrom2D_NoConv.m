function Cell=DoGfrom2D_NoConv(Em)
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

ConesToCenter=floor(0.3736*Em^2 + 0.01126*Em + 1.026);
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
Lc=LcAmp*(normpdf(x,0,LcSigma));

McSigma=CenterSigma;
McAmp=crM;
Mc=McAmp*(normpdf(x,0,McSigma));

ScSigma=0;
ScAmp=crS;
Sc=ScAmp*(normpdf(x,0,ScSigma));

LsSigma=SurroundSigma;
LsAmp=srL;
Ls=LsAmp*(normpdf(x,0,LsSigma));

MsSigma=SurroundSigma;
MsAmp=srM;
Ms=MsAmp*(normpdf(x,0,MsSigma));

SsSigma=0;
SsAmp=crS;
Ss=SsAmp*(normpdf(x,0,SsSigma));

% Total cell response:
CellResp=Lc+Ls+Mc+Ms+Sc+Ss;

%% FFT OF THE RESPONSE PROFILE

TotalLResp=Lc+Ls;
TotalMResp=Mc+Ms;
LMSumResp=TotalLResp+TotalMResp;
LMDiffResp=TotalLResp-TotalMResp;

TotalLFFT=fftshift(fft(TotalLResp,w));  %*1/Fs
TotalMFFT=fftshift(fft(TotalMResp,w));
LMSumFFT=fftshift(fft(LMSumResp,w));
LMDiffFFT=fftshift(fft(LMDiffResp,w));

TotalLPower=abs(TotalLFFT);
TotalMPower=abs(TotalMFFT);
LMSumPower=abs(LMSumFFT);
LMDiffPower=abs(LMDiffFFT);

TotalLphase=angle(TotalLFFT)*180/pi;
TotalMphase=angle(TotalMFFT)*180/pi;
LMSumphase=angle(LMSumFFT)*180/pi;
LMDiffphase=angle(LMDiffFFT)*180/pi;

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
    
    %FFT sampling indices
    X1=SpatFreq;
    X2=-SpatFreq;

    FindIndex=[X1 X2];
    for b=1:length(FindIndex)
        if b>0
            [~, minInd(b)] = min(abs(Freq-FindIndex(b)));
        elseif b<0
            [~, minInd(b)] = min(abs(Freq-FindIndex(b)));
        end
    end
    
    %%%%%%% Get FFT values

    lAmp=TotalLPower(minInd(1));
    mAmp=TotalMPower(minInd(1));
    sumAmp=LMSumPower(minInd(1));
    diffAmp=LMDiffPower(minInd(1));
    normAmp=(max(TotalLPower)+max(TotalMPower))/100;
    
    lPhase=TotalLphase(minInd(1));
    mPhase=TotalMphase(minInd(1));
    sumPhase=LMSumphase(minInd(1));
    diffPhase=LMDiffphase(minInd(1));
    
    TotalLData.Amplitude(t)=lAmp/normAmp;
    TotalMData.Amplitude(t)=mAmp/normAmp;
    LMSumData.Amplitude(t)=sumAmp/normAmp;
    LMDiffData.Amplitude(t)=diffAmp/normAmp;
    
    TotalLData.Phase(t)=lPhase;
    TotalMData.Phase(t)=mPhase;
    LMSumData.Phase(t)=sumPhase;
    LMDiffData.Phase(t)=diffPhase;

end

%% TEST FOR CHROMATIC OPPONENCY 

%Assess chromatic opponency using the phase relationship
AmpIndex=LMSumData.Amplitude(1)-LMDiffData.Amplitude(1);
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
Cell.Eccentricity=Em;

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
Cell.ResponseFunctions.LMSumResponse.Amplitude=LMSumData.Amplitude;
Cell.ResponseFunctions.LMSumResponse.Phase=LMSumData.Phase;
Cell.ResponseFunctions.LMDiffResponse.Amplitude=LMDiffData.Amplitude;
Cell.ResponseFunctions.LMDiffResponse.Phase=LMDiffData.Phase;

%% PLOTTING

%PlotASingleCell;

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
% Ymax=max([LMSumData.Amplitude]);
% if Ymax==0
%     Ymax=10;
% end
% PlotTotAmp=loglog(Freqs,LMSumData.Amplitude,'ro');
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
% PlotTotPhase=semilogx(Freqs,abs(LMSumData.Phase),'ro');
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
