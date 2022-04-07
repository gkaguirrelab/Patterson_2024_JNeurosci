function Cell=Grid2DDoG_CrysS_H2(Em,Universe,thetas,SpatFreqs,SFLines,OriCircles)
%% GENERATE CELL
Cell=GenerateCell_CrysS_H2(Em);

%% MODEL UNIVERSE

NumDegs=64; %Number of degrees across a dimension (to create an NumDeg-by-NumDeg patch of retina)
XRes=2^13; %How resolved is the patch? (i.e., matrix length = XRes)
Fs=XRes/NumDegs; %Sampling frequency
X = -NumDegs/2:1/Fs:NumDegs/2-(1/Fs); %Where in space is this NumDeg-by-NumDeg patch?
W=length(X); %Determine length of X dimension for FFT frequency determination
Freq=(-W/2:(W/2)-1)*(Fs/W); %Frequency domain values for FFT

Universe.NumDegs=NumDegs;
Universe.Fs=Fs;
Universe.X=X;
Universe.W=W;
Universe.Freq=Freq;

%% GENERATE GRATING FFT EVALUATION LINES 
%This only has to be done once and the indices can be used to fetch
%amp/phase for all cell RF FFTs
thetas=[0 30 60 90 120 150];
% thetas=[0];
SpatFreqs=[1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32];

SFLines=ComputeSFLine(thetas,W,Freq);
OriCircles=ComputeOriCircle(SpatFreqs,W,Freq);

%% IMPORT DATA & CONVERT TO DEGREES

LCenterCoords=Cell.Center.LCones.Coords;
LCenterCoordsDeg=(LCenterCoords*1000)/200;
LCenterCoordsDeg(isnan(LCenterCoordsDeg(:,1)),:) = [];
LCenterKernel=(Cell.Center.LCones.EachWeight);
LCenterKernel(LCenterKernel==0)=[];

MCenterCoords=Cell.Center.MCones.Coords;
MCenterCoordsDeg=(MCenterCoords*1000)/200;
MCenterCoordsDeg(isnan(MCenterCoordsDeg(:,1)),:) = [];
MCenterKernel=(Cell.Center.MCones.EachWeight);
MCenterKernel(MCenterKernel==0)=[];

SCenterCoords=Cell.Center.SCones.Coords;
SCenterCoordsDeg=(SCenterCoords*1000)/200;
SCenterCoordsDeg(isnan(SCenterCoordsDeg(:,1)),:) = [];
SCenterKernel=(Cell.Center.SCones.EachWeight);
SCenterKernel(SCenterKernel==0)=[];

LSurroundCoords=Cell.Surround.LCones.Coords;
LSurroundCoordsDeg=(LSurroundCoords*1000)/200;
LSurroundCoordsDeg(isnan(LSurroundCoordsDeg(:,1)),:) = [];
LSurroundKernel=(Cell.Surround.LCones.EachWeight);
LSurroundKernel(LSurroundKernel==0)=[];

MSurroundCoords=Cell.Surround.MCones.Coords;
MSurroundCoordsDeg=(MSurroundCoords*1000)/200;
MSurroundCoordsDeg(isnan(MSurroundCoordsDeg(:,1)),:) = [];
MSurroundKernel=(Cell.Surround.MCones.EachWeight);
MSurroundKernel(MSurroundKernel==0)=[];

SSurroundCoords=Cell.Surround.SCones.Coords;
SSurroundCoordsDeg=(SSurroundCoords*1000)/200;
SSurroundCoordsDeg(isnan(SSurroundCoordsDeg(:,1)),:) = [];
SSurroundKernel=(Cell.Surround.SCones.EachWeight);
SSurroundKernel(SSurroundKernel==0)=[];

CenterRadius=max(Cell.Center.CenterConeDistances);
CenterRadiusDeg=(CenterRadius*1000)/200;

SurroundRadius=max(Cell.Surround.SurroundConeDistances);
SurroundRadiusDeg=(SurroundRadius*1000)/200;

%% MODEL UNIVERSE

%Unpack the universe variables
v2struct(Universe);

%% RESPONSE PROFILE GENERATION

%Assign a sigma to each cone based on photoreceptor diameter
ConeRad=3.955*exp(0.0163*Em)-3.149*exp(-1.288*Em); %From Fig 8, Packer et al, 1989 JCN. THIS IS IN MICRONS
ConeRadDeg=(ConeRad)/200;
ConeSigma=ConeRadDeg;

%How strong is the surround with respect to the center?
CSStrengthRatio=Cell.CSStrengthRatio;

%Define a truncated RF space to generate the Gaussians. These will be
%padded later
PaddingRadius=2.5;
RFX=-PaddingRadius:1/Fs:PaddingRadius-1/Fs;
[RFXm, RFYm] = meshgrid(RFX, RFX);

%CENTER RESPONSE PROFILES

%Apply a Gaussian to each center L cone coordinate...
if size(LCenterCoordsDeg,1)>=1
    for t=1:size(LCenterCoordsDeg,1)
    CgaussL(:,:,t)=LCenterKernel(t)*exp(-((RFXm-LCenterCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-LCenterCoordsDeg(t,2)).^2/2/ConeSigma^2));
    end
    AllLc=sum(CgaussL,3);
else
    AllLc=zeros(length(RFXm));
end

%...and for each center M cone coordinate
if size(MCenterCoordsDeg,1)>=1
    for t=1:size(MCenterCoordsDeg,1)
    CgaussM(:,:,t)=MCenterKernel(t)*exp(-((RFXm-MCenterCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-MCenterCoordsDeg(t,2)).^2/2/ConeSigma^2));
    end
    AllMc=sum(CgaussM,3);
else
    AllMc=zeros(length(RFXm));
end

%...and for each center S cone coordinate
if size(SCenterCoordsDeg,1)>=1
    for t=1:size(SCenterCoordsDeg,1)
    CgaussS(:,:,t)=SCenterKernel(t)*exp(-((RFXm-SCenterCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-SCenterCoordsDeg(t,2)).^2/2/ConeSigma^2));
    end
    AllSc=sum(CgaussS,3);
else
    AllSc=zeros(length(RFXm));
end

%SURROUND RESPONSE PROFILES

%Apply a Gaussian to each surround L cone coordinate...
for t=1:length(LSurroundCoordsDeg)
SgaussL(:,:,t)=LSurroundKernel(t)*exp(-((RFXm-LSurroundCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-LSurroundCoordsDeg(t,2)).^2/2/ConeSigma^2));
end
AllLs=-sum(SgaussL,3)*CSStrengthRatio;

%...and for each surround M cone coordinate
for t=1:length(MSurroundCoordsDeg)
SgaussM(:,:,t)=MSurroundKernel(t)*exp(-((RFXm-MSurroundCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-MSurroundCoordsDeg(t,2)).^2/2/ConeSigma^2));
end
AllMs=-sum(SgaussM,3)*CSStrengthRatio;

%...and for each surround S cone coordinate
for t=1:length(SSurroundCoordsDeg)
SgaussS(:,:,t)=SSurroundKernel(t)*exp(-((RFXm-SSurroundCoordsDeg(t,1)).^2/2/ConeSigma^2 + (RFYm-SSurroundCoordsDeg(t,2)).^2/2/ConeSigma^2));
end
AllSs=-sum(SgaussS,3)*CSStrengthRatio;

%The computed L/L+M proportions should be identical to what
%'GenerateCell.m' outputs. 
cpM=(sum(sum(AllMc)))/(sum(sum(AllLc))+sum(sum(AllMc))+sum(sum(AllSc)));
cpL=(sum(sum(AllLc)))/(sum(sum(AllLc))+sum(sum(AllMc))+sum(sum(AllSc)));
cpS=(sum(sum(AllSc)))/(sum(sum(AllLc))+sum(sum(AllMc))+sum(sum(AllSc)));

spM=(sum(sum(AllMs)))/(sum(sum(AllLs))+sum(sum(AllMs))+sum(sum(AllSs)));
spL=(sum(sum(AllLs)))/(sum(sum(AllLs))+sum(sum(AllMs))+sum(sum(AllSs)));
spS=(sum(sum(AllSs)))/(sum(sum(AllLs))+sum(sum(AllMs))+sum(sum(AllSs)));

%Form the center-surround DoG
AllL=AllLc+AllLs;
AllM=AllMc+AllMs;
AllS=AllSc+AllSs;

%Now pad the arrays to make them the same size as the patch
AllLc=padarray(AllLc,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllMc=padarray(AllMc,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllSc=padarray(AllSc,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);

AllLs=padarray(AllLs,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllMs=padarray(AllMs,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllSs=padarray(AllSs,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);

AllL=padarray(AllL,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllM=padarray(AllM,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllS=padarray(AllS,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);

%% FFT

TotalLResp=AllLc+AllLs;
TotalMResp=AllMc+AllMs;
TotalSResp=AllSc+AllSs;

LMSumResp=TotalLResp+TotalMResp;
LMDiffResp=TotalLResp-TotalMResp;

% LcFFT=fftshift(fft2(ifftshift(AllLc)));
% Lcpower=abs(LcFFT);
% Lcphase=angle(LcFFT)*180/pi;
% 
% LsFFT=fftshift(fft2(ifftshift(AllLs)));
% Lspower=abs(LsFFT);
% Lsphase=angle(LsFFT)*180/pi;
% 
% McFFT=fftshift(fft2(ifftshift(AllMc)));
% Mcpower=abs(McFFT);
% Mcphase=angle(McFFT)*180/pi;
% 
% MsFFT=fftshift(fft2(ifftshift(AllMs)));
% Mspower=abs(MsFFT);
% Msphase=angle(MsFFT)*180/pi;

TotalLFFT=fftshift(fft2(ifftshift(TotalLResp)));
TotalLpower=abs(TotalLFFT);
TotalLphase=angle(TotalLFFT)*180/pi;

TotalMFFT=fftshift(fft2(ifftshift(TotalMResp)));
TotalMpower=abs(TotalMFFT);
TotalMphase=angle(TotalMFFT)*180/pi;

TotalSFFT=fftshift(fft2(ifftshift(TotalSResp)));
TotalSpower=abs(TotalSFFT);
TotalSphase=angle(TotalSFFT)*180/pi;

LMSumFFT=fftshift(fft2(ifftshift(LMSumResp)));
LMSumpower=abs(LMSumFFT);
LMSumphase=angle(LMSumFFT)*180/pi;

LMDiffFFT=fftshift(fft2(ifftshift(LMDiffResp)));
LMDiffpower=abs(LMDiffFFT);
LMDiffphase=angle(LMDiffFFT)*180/pi;

%% EVALUATE FFT RESPONSES GIVEN SF LINES OR ORI CIRCLES
for t=1:length(thetas)
    OrientationLabel=strcat('Theta',num2str(thetas(t)),'deg');
    %Read the appropriate values from the RF FFTs
    lAmp=TotalLpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    mAmp=TotalMpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    sAmp=TotalSpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    sumAmp=LMSumpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    diffAmp=LMDiffpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    normAmp=(max(max(TotalLpower))+max(max(TotalMpower))+max(max(TotalSpower)))/100;

    lPhase=TotalLphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    mPhase=TotalMphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    sPhase=TotalSphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    sumPhase=LMSumphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    diffPhase=LMDiffphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);

    TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=lAmp/normAmp;
    TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=mAmp/normAmp;
    TotalSData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=sAmp/normAmp;
    LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=sumAmp/normAmp;
    LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=diffAmp/normAmp;

    TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Phase=lPhase;%unwrap(lPhase);
    TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Phase=mPhase;%unwrap(mPhase);
    TotalSData.(matlab.lang.makeValidName(OrientationLabel)).Phase=sPhase;%unwrap(sPhase);
    LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Phase=sumPhase;%unwrap(sumPhase);
    LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Phase=diffPhase;%unwrap(diffPhase);
end

for s=1:length(SpatFreqs)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    lAmp=TotalLpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mAmp=TotalMpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sAmp=TotalSpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sumAmp=LMSumpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    diffAmp=LMDiffpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    normAmp=(max(max(TotalLpower))+max(max(TotalMpower))+max(max(TotalSpower)))/100;

    lPhase=TotalLphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mPhase=TotalMphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sPhase=TotalSphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sumPhase=LMSumphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    diffPhase=LMDiffphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);

    TotalLData.(matlab.lang.makeValidName(SFLabel)).Amplitude=lAmp/normAmp;
    TotalMData.(matlab.lang.makeValidName(SFLabel)).Amplitude=mAmp/normAmp;
    TotalSData.(matlab.lang.makeValidName(SFLabel)).Amplitude=sAmp/normAmp;
    LMSumData.(matlab.lang.makeValidName(SFLabel)).Amplitude=sumAmp/normAmp;
    LMDiffData.(matlab.lang.makeValidName(SFLabel)).Amplitude=diffAmp/normAmp;

    TotalLData.(matlab.lang.makeValidName(SFLabel)).Phase=lPhase;%unwrap(lPhase);
    TotalMData.(matlab.lang.makeValidName(SFLabel)).Phase=mPhase;%unwrap(mPhase);
    TotalSData.(matlab.lang.makeValidName(SFLabel)).Phase=sPhase;%unwrap(sPhase);
    LMSumData.(matlab.lang.makeValidName(SFLabel)).Phase=sumPhase;%unwrap(sumPhase);
    LMDiffData.(matlab.lang.makeValidName(SFLabel)).Phase=diffPhase;%unwrap(diffPhase);
end

%% TEST FOR CHROMATIC OPPONENCY 

%Assess chromatic opponency using the phase relationship
AmpIndex=LMSumData.Theta0deg.Amplitude(1)-LMDiffData.Theta0deg.Amplitude(1);
PhaseIndex=abs(TotalLData.Theta0deg.Phase(1))-abs(TotalMData.Theta0deg.Phase(1));
sPhaseIndex=abs(TotalSData.Theta0deg.Phase(1));

if AmpIndex<0 %L-M response is greater than L+M response at full-field stimulation
    if PhaseIndex<-175 %L phase = 0, M phase = 180
        ChromTag='L-dominated';
    elseif PhaseIndex>175 %M phase = 0, L phase = 180
        ChromTag='M-dominated';
    end
elseif AmpIndex>0 %L+M response is greater than/equal to L-M response at full-field stimulation
    if abs(sPhaseIndex-PhaseIndex)>175
        ChromTag='-S vs +LM';
    else
        ChromTag='Achromatic';
    end
end

%% SAVE DATA

Cell.ChromTag=ChromTag;
Cell.Eccentricity=Em;

Cell.RFProfile.DFRadiusDeg=CenterRadiusDeg;
Cell.RFProfile.SpatialFrequencies=SFLines;

% Cell.RFProfile.LResponse.Center=AllLc;
Cell.RFProfile.LResponse.CenterAmplitude=cpL;
% Cell.RFProfile.LResponse.Surround=AllLs;
Cell.RFProfile.LResponse.SurroundAmplitude=spL;

% Cell.RFProfile.MResponse.Center=AllMc;
Cell.RFProfile.MResponse.CenterAmplitude=cpM;
% Cell.RFProfile.MResponse.Surround=AllMs;
Cell.RFProfile.MResponse.SurroundAmplitude=spM;

% Cell.RFProfile.SResponse.Center=AllSc;
Cell.RFProfile.SResponse.CenterAmplitude=cpS;
% Cell.RFProfile.SResponse.Surround=AllSs;
Cell.RFProfile.SResponse.SurroundAmplitude=spS;

Cell.ResponseFunctions.LIsoResponse=TotalLData;
Cell.ResponseFunctions.MIsoResponse=TotalMData;
Cell.ResponseFunctions.SIsoResponse=TotalSData;
Cell.ResponseFunctions.LMSumResponse=LMSumData;
Cell.ResponseFunctions.LMDiffResponse=LMDiffData;

%% PLOT THE RFs

% %%%%%%%%%%%% RF and FFT
% 
% set(0,'DefaultFigureWindowStyle','docked');
% FiguresOpen=length(findobj(0,'type','figure'));
% figure(FiguresOpen+1);
% 
% subplot(2,2,1);
% imagesc(X,X,AllLc+AllLs);
% title('L-Cone Response Profile');
% xlabel('x (deg)');
% ylabel('y (deg)'); 
% axis square;
% axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);
% 
% subplot(2,2,2);
% imagesc(X,X,AllMc+AllMs);
% title('M-Cone Response Profile');
% xlabel('x (deg)');
% ylabel('y (deg)'); 
% axis square;
% axis([-SurroundRadiusDeg*2 SurroundRadiusDeg*2 -SurroundRadiusDeg*2 SurroundRadiusDeg*2]);
% 
% subplot(2,2,3);
% imagesc(Freq,Freq,TotalLpower);
% title('L-Cone Frequency Domain');
% xlabel('Frequency (cpd)');
% ylabel('Frequency (cpd)'); 
% axis square;
% 
% subplot(2,2,4);
% imagesc(Freq,Freq,TotalMpower);
% title('M-Cone Frequency Domain');
% xlabel('Frequency (cpd)');
% ylabel('Frequency (cpd)'); 
% axis square;
% 
%     %%%%%%%%%%%%% DoG OUTPUTS
%     
% for f=1:length(thetas)
%     OrientationLabel=strcat('Theta',num2str(thetas(f)),'deg');
%     Freqs=Cell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).zFreqs;
%     Line=Cell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).zIndices;
%     PlotL=TotalLData.(matlab.lang.makeValidName(OrientationLabel));
%     PlotM=TotalMData.(matlab.lang.makeValidName(OrientationLabel));
%     PlotLMSum=LMSumData.(matlab.lang.makeValidName(OrientationLabel));
%     PlotLMDiff=LMDiffData.(matlab.lang.makeValidName(OrientationLabel));
%     
%     figure(f+1+FiguresOpen);
%     centerstrength=strcat('Lc=',num2str(cpL*100),', Mc=',num2str(cpM*100));
%     surroundstrength=strcat('Ls=',num2str(spL*100),', Ms=',num2str(spM*100));
%     ori=strcat('Theta',num2str(thetas(f)),'deg');
% 
%     subplot(2,3,1);
%     Ymax=max([PlotL.Amplitude PlotM.Amplitude]);
%     Ymin=min([PlotL.Amplitude PlotM.Amplitude]);
%     if Ymax==0
%         Ymax=10;
%     end
%     PlotLAmp=loglog(Freqs,PlotL.Amplitude,'ro');
%     hold on;
%     PlotMAmp=loglog(Freqs,PlotM.Amplitude,'go');
%     %PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Response amplitude (a.u.)'); 
%     axis([0.005 15 .1 Ymax*10]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     legend([PlotLAmp,PlotMAmp],'Location','southwest')
%     set(PlotLAmp,...
%         'DisplayName','L',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[204 0 0]/255)
%     set(PlotMAmp,...
%         'DisplayName','M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[119 172 48]/255)
%     %     set(PlotSAmp,...
%     %         'DisplayName','S',...
%     %         'LineWidth',.5,...
%     %         'LineStyle',':',...
%     %         'Color','k',...
%     %         'MarkerFaceColor',[0 0.45 0.74])
% 
%     subplot(2,3,4);
%     PlotLPhase=semilogx(Freqs,abs(PlotL.Phase),'ro');
%     hold on;
%     PlotMPhase=semilogx(Freqs,abs(PlotM.Phase),'go');
%     %PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Phase (degrees)'); 
%     axis([0.005 15 -20 200]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     %legend([PlotLPhase,PlotMPhase],'Location','northeast')
%     set(PlotLPhase,...
%         'DisplayName','L',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[204 0 0]/255)
%     set(PlotMPhase,...
%         'DisplayName','M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[119 172 48]/255)
%     %     set(PlotSPhase,...
%     %         'DisplayName','S',...
%     %         'LineWidth',.5,...
%     %         'LineStyle',':',...
%     %         'Color','k',...
%     %         'MarkerFaceColor',[0 0.45 0.74])
% 
%     subplot(2,3,2);
%     Ymax=max([PlotLMSum.Amplitude]);
%     Ymin=min([PlotLMSum.Amplitude]);
%     if Ymax==0
%         Ymax=10;
%     end
%     PlotTotAmp=loglog(Freqs,PlotLMSum.Amplitude,'ro');
%     hold on;
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Response amplitude (a.u.)'); 
%     axis([0.005 15 .1 Ymax*10]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     legend(PlotTotAmp,'Location','southwest')
%     set(PlotTotAmp,...
%         'DisplayName','L+M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[.8 .8 .8])
% 
%     subplot(2,3,5);
%     PlotTotPhase=semilogx(Freqs,abs(PlotLMSum.Phase),'ro');
%     hold on;
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Phase (degrees)'); 
%     axis([0.005 15 -20 200]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     %legend(PlotTotPhase,'Location','northeast')
%     set(PlotTotPhase,...
%         'DisplayName','L+M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[.8 .8 .8])
% 
%     subplot(2,3,3);
%     Ymax=max([PlotLMDiff.Amplitude]);
%     Ymin=min([PlotLMDiff.Amplitude]);
%     if Ymax==0
%         Ymax=10;
%     end
%     PlotLMDiffAmp=loglog(Freqs,PlotLMDiff.Amplitude,'ro');
%     hold on;
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Response amplitude (a.u.)'); 
%     axis([0.005 15 .1 Ymax*10]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     legend(PlotLMDiffAmp,'Location','southwest')
%     set(PlotLMDiffAmp,...
%         'DisplayName','L-M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[.31 .31 .31])
% 
%     subplot(2,3,6);
%     PlotLMDiffPhase=semilogx(Freqs,abs(PlotLMDiff.Phase),'ro');
%     hold on;
%     %title({centerstrength;surroundstrength});
%     xlabel('Stimulus frequency (cpd)');
%     ylabel('Phase (degrees)'); 
%     axis([0.005 15 -20 200]);
%     set(gca,'TickDir','in','TickLength', [.005 .005]);box off
%     %legend(PlotTotPhase,'Location','northeast')
%     set(PlotLMDiffPhase,...
%         'DisplayName','L-M',...
%         'LineWidth',.5,...
%         'LineStyle',':',...
%         'Color','k',...
%         'MarkerFaceColor',[.31 .31 .31])
% 
%     axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     text(0.5, 1,{centerstrength;surroundstrength;ori},'HorizontalAlignment','center','VerticalAlignment', 'top');
% end
% set(0,'DefaultFigureWindowStyle','default');

%% WORKSHOP
% MAKE A STIMULUS GRATING

% Freqs=[1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32];
% thetas=[0 30 60 90 120 150];
% SFLabel=strcat('SF=',num2str(Freqs(t)));
% t=2;
% for t=1:length(Freqs)
%     for th=1:length(thetas)
%         OrientationLabel=strcat('Theta',num2str(thetas(th)),'deg');
%         SpatFreq=Freqs(t);
%         theta = thetas(3); %Grating orientation (degrees)
%         phase = 0; %Phase (degree)
%         phi = (phase/360)*2*pi;
%         thetaRad=(theta/360)*2*pi; % convert theta (orientation) to radians
%         Xt = Xm*cos(thetaRad); %Compute proportion of Xm for given orientation
%         Yt = Ym*sin(thetaRad); %Compute proportion of Ym for given orientation
%         XYt=[Xt+Yt]; %Sum X and Y components
%         grating = cos(2*pi*SpatFreq*XYt+phi);
% 
%         %Plot it if you want
%         figure;imagesc(X,X,grating);axis square;
% 
%         %FFT sampling indices
%         X1=SpatFreq*cos(thetaRad);
%         Y1=SpatFreq*sin(thetaRad);
%         X2=-SpatFreq*cos(thetaRad);
%         Y2=-SpatFreq*sin(thetaRad);
% 
%         FindIndex(th,:)=[X1 Y1 X2 Y2];
%         for b=1:length(FindIndex)
%             if b>0
%                 [~, minInd(th,b)] = min(abs(Freq-FindIndex(b)));
%             elseif b<0
%                 [~, minInd(th,b)] = min(abs(Freq-FindIndex(b)));
%             end
%         end
% 
%         %Read the appropriate values from the RF FFTs
%         lAmp=TotalLpower(minInd(th,1),minInd(th,2));
%         mAmp=TotalMpower(minInd(th,1),minInd(th,2));
%         sumAmp=LMSumpower(minInd(th,1),minInd(th,2));
%         diffAmp=LMDiffpower(minInd(th,1),minInd(th,2));
%         normAmp=(max(max(TotalLpower))+max(max(TotalMpower)))/100;
%         
%         lPhase=TotalLphase(minInd(th,1),minInd(th,2));
%         mPhase=TotalMphase(minInd(th,1),minInd(th,2));
%         sumPhase=LMSumphase(minInd(th,1),minInd(th,2));
%         diffPhase=LMDiffphase(minInd(th,1),minInd(th,2));
%         
%         TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(t)=lAmp/normAmp;
%         TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(t)=mAmp/normAmp;
%         LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(t)=sumAmp/normAmp;
%         LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude(t)=diffAmp/normAmp;
%         
%         TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Phase(t)=lPhase;
%         TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Phase(t)=mPhase;
%         LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Phase(t)=sumPhase;
%         LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Phase(t)=diffPhase;
%     end
% 
% end

% ORI CIRCLES

% SpatFreq=2;
% 
% [~, radInd]=min(abs(Freq-SpatFreq)); %Find the index value corresponding to SpatFreq
% [~, oriInd]=min(abs(Freq)); %Find the index value corresponding to SpatFreq=0 (i.e., origin)
% 
% CircleSubstrate=zeros(W);
% Radius=abs(radInd-oriInd);
% 
% circle=MidpointCircle(CircleSubstrate,Radius,oriInd,oriInd,1);
% imagesc(circle);
% 
% SFLcirc=TotalLpower.*circle;
% [Li,Lj,Lamp] = find(SFLcirc);
% [~,LmaxInd]=max(Lamp);
% LthetaInd=[Li(LmaxInd),Lj(LmaxInd)];
% SFMcirc=TotalMpower.*circle;
% [Mi,Mj,Mamp] = find(SFMcirc);
% [~,MmaxInd]=max(Mamp);
% 
% LiFreq=Freq(LthetaInd(1));
% LjFreq=Freq(LthetaInd(2));
% 
% PrefOriFreq=atand(LiFreq/LjFreq);
% if PrefOriFreq<0
%     PrefOriFreq=180+atand(LiFreq/LjFreq);
% end
% 
% %%
% %"Logsample" the dataset for plotting
% for z=1:12
%     zInd(z,1)=2^(z-1);
% end
% 
% figure;
% loglog(Freqs,TotalLpower(Line(zInd)),'r.');
% hold on
% loglog(Freqs,TotalMpower(Line(zInd)),'g.');
% 
% figure;
% semilogx(Freqs,TotalLphase(Line(zInd)),'r.');
% hold on
% semilogx(Freqs,TotalMphase(Line(zInd)),'g.');

