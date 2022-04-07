%% IMPORT DATA STRUCTURE

% Data=load('/Volumes/Macintosh HDD/Documents/SCHOOL/RF Modeling Data/2DMidgetRun_24-Nov-2015/2DMidgetRun_24-Nov-2015.mat');
% CellList=fieldnames(Data);

% Identifying information
% ImportCell=Data.Ecc1mm.Cell5;
% CellName='Data.Ecc1mm.Cell5';
ImportCell=Data.(matlab.lang.makeValidName(Cells{2810}));
%CellName='Data.Ecc7_2322mm';

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

%% LOAD VARIABLES

Em=ImportCell.Eccentricity;

cpL=ImportCell.CenterWeights(1);
cpM=ImportCell.CenterWeights(2);
cpS=ImportCell.CenterWeights(3);

spL=ImportCell.SurroundWeights(1);
spM=ImportCell.SurroundWeights(2);
spS=ImportCell.SurroundWeights(3);

LCenterCoords=ImportCell.Center.LCones.Coords;
LCenterCoordsDeg=(LCenterCoords*1000)/200;
LCenterCoordsDeg(isnan(LCenterCoordsDeg(:,1)),:) = [];
LCenterKernel=(~isnan(ImportCell.Center.LCones.Coords(:,1))).*(ImportCell.Center.CenterKernel);
LCenterKernel(LCenterKernel==0)=[];

MCenterCoords=ImportCell.Center.MCones.Coords;
MCenterCoordsDeg=(MCenterCoords*1000)/200;
MCenterCoordsDeg(isnan(MCenterCoordsDeg(:,1)),:) = [];
MCenterKernel=(~isnan(ImportCell.Center.MCones.Coords(:,1))).*(ImportCell.Center.CenterKernel);
MCenterKernel(MCenterKernel==0)=[];

LSurroundCoords=ImportCell.Surround.LCones.Coords;
LSurroundCoordsDeg=(LSurroundCoords*1000)/200;
LSurroundCoordsDeg(isnan(LSurroundCoordsDeg(:,1)),:) = [];
LSurroundKernel=(~isnan(ImportCell.Surround.LCones.Coords(:,1))).*(ImportCell.Surround.SurroundKernel);
LSurroundKernel(LSurroundKernel==0)=[];

MSurroundCoords=ImportCell.Surround.MCones.Coords;
MSurroundCoordsDeg=(MSurroundCoords*1000)/200;
MSurroundCoordsDeg(isnan(MSurroundCoordsDeg(:,1)),:) = [];
MSurroundKernel=(~isnan(ImportCell.Surround.MCones.Coords(:,1))).*(ImportCell.Surround.SurroundKernel);
MSurroundKernel(MSurroundKernel==0)=[];

CenterRadius=max(ImportCell.Center.CenterConeDistances);
CenterRadiusDeg=(CenterRadius*1000)/200;

SurroundRadius=max(ImportCell.Surround.SurroundConeDistances);
SurroundRadiusDeg=(SurroundRadius*1000)/200;

TotalLData=ImportCell.ResponseFunctions.LIsoResponse;
TotalMData=ImportCell.ResponseFunctions.MIsoResponse;
LMSumData=ImportCell.ResponseFunctions.LMSumResponse;
LMDiffData=ImportCell.ResponseFunctions.LMDiffResponse;

% %CD to plot directory
% PlotDirectory=strcat('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid/Analysis/Cell Plots/',CellName);
% if (exist(PlotDirectory) == 0)
%     mkdir(PlotDirectory);
%     cd(PlotDirectory);
% else
%     cd(PlotDirectory);
% end

%% RESPONSE PROFILE GENERATION

%Assign a sigma to each cone based on photoreceptor diameter
ConeRad=3.955*exp(0.0163*Em)-3.149*exp(-1.288*Em); %From Fig 8, Packer et al, 1989 JCN. THIS IS IN MICRONS
ConeRadDeg=(ConeRad)/200;
ConeSigma=ConeRadDeg;

%How strong is the surround with respect to the center?
CSStrengthRatio=ImportCell.CSStrengthRatio;

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

%The computed L/L+M proportions should be identical to what
%'GenerateCell.m' outputs. 
cpM=(sum(sum(AllMc)))/(sum(sum(AllLc))+sum(sum(AllMc)));
cpL=(sum(sum(AllLc)))/(sum(sum(AllLc))+sum(sum(AllMc)));

spM=(sum(sum(AllMs)))/(sum(sum(AllLs))+sum(sum(AllMs)));
spL=(sum(sum(AllLs)))/(sum(sum(AllLs))+sum(sum(AllMs)));

%Form the center-surround DoG
AllL=AllLc+AllLs;
AllM=AllMc+AllMs;

%Now pad the arrays to make them the same size as the patch
AllLc=padarray(AllLc,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllMc=padarray(AllMc,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllLs=padarray(AllLs,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllMs=padarray(AllMs,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllL=padarray(AllL,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);
AllM=padarray(AllM,[(NumDegs/2-PaddingRadius)*Fs (NumDegs/2-PaddingRadius)*Fs]);

%% FFT

TotalLResp=AllLc+AllLs;
TotalMResp=AllMc+AllMs;
LMSumResp=TotalLResp+TotalMResp;
LMDiffResp=TotalLResp-TotalMResp;

TotalLFFT=fftshift(fft2(ifftshift(TotalLResp)));
TotalLpower=abs(TotalLFFT);
TotalLphase=angle(TotalLFFT)*180/pi;

TotalMFFT=fftshift(fft2(ifftshift(TotalMResp)));
TotalMpower=abs(TotalMFFT);
TotalMphase=angle(TotalMFFT)*180/pi;

%% PLOT THE RFs

%%%%%%%%%%%% RF and FFT

set(0,'DefaultFigureWindowStyle','docked');
FiguresOpen=length(findobj(0,'type','figure'));
figure(FiguresOpen+1);

subplot(2,2,1);
imagesc(X,-X,TotalLResp);
title('L-Cone Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*1.1 SurroundRadiusDeg*1.1 -SurroundRadiusDeg*1.1 SurroundRadiusDeg*1.1]);

subplot(2,2,2);
imagesc(X,-X,TotalMResp);
title('M-Cone Response Profile');
xlabel('x (deg)');
ylabel('y (deg)'); 
axis square;
axis([-SurroundRadiusDeg*1.1 SurroundRadiusDeg*1.1 -SurroundRadiusDeg*1.1 SurroundRadiusDeg*1.1]);

subplot(2,2,3);
imagesc(Freq,-Freq,TotalLpower);
title('L-Cone Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;
axis([-5 5 -5 5]);

subplot(2,2,4);
imagesc(Freq,-Freq,TotalMpower);
title('M-Cone Frequency Domain');
xlabel('Frequency (cpd)');
ylabel('Frequency (cpd)'); 
axis square;
axis([-5 5 -5 5]);

% colormap jet;

%% %%%%%%%%%%% DoG OUTPUTS

set(0,'DefaultFigureWindowStyle','docked');
FiguresOpen=length(findobj(0,'type','figure'));

thetas=[0 30 60 90 120 150];   
    
for f=1:length(thetas)
    OrientationLabel=strcat('Theta',num2str(thetas(f)),'deg');
    Freqs=ImportCell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).Freqs;
    Line=ImportCell.RFProfile.SpatialFrequencies.(matlab.lang.makeValidName(OrientationLabel)).Indices;
    PlotL=TotalLData.(matlab.lang.makeValidName(OrientationLabel));
    PlotM=TotalMData.(matlab.lang.makeValidName(OrientationLabel));
%     PlotS=TotalSData.(matlab.lang.makeValidName(OrientationLabel));
    PlotLMSum=LMSumData.(matlab.lang.makeValidName(OrientationLabel));
    PlotLMDiff=LMDiffData.(matlab.lang.makeValidName(OrientationLabel));
    
    Ymax=max([PlotL.Amplitude PlotM.Amplitude PlotLMSum.Amplitude PlotLMDiff.Amplitude]);
    Ymin=min([PlotL.Amplitude PlotM.Amplitude PlotLMSum.Amplitude PlotLMDiff.Amplitude]);
    
    figure(f+FiguresOpen);
    %set(gcf,'Position',[1 1 845 325]);
    centerstrength=strcat('Lc=',num2str(cpL*100),', Mc=',num2str(cpM*100));
    surroundstrength=strcat('Ls=',num2str(spL*100),', Ms=',num2str(spM*100));
    ori=strcat('Theta',num2str(thetas(f)),'deg');

    subplot(4,3,[1 7]);
%     Ymax=max([PlotL.Amplitude PlotM.Amplitude]);
%     Ymin=min([PlotL.Amplitude PlotM.Amplitude]);
    if Ymax==0
        Ymax=10;
    end
    PlotLAmp=loglog(Freqs,PlotL.Amplitude,'r');
    hold on;
    PlotMAmp=loglog(Freqs,PlotM.Amplitude,'g');
%     PlotSAmp=loglog(Freqs,PlotS.Amplitude,'b');
    %title({centerstrength;surroundstrength});
    %xlabel('Stimulus frequency (cpd)');
    ylabel('Response amplitude (a.u.)'); 
    axis([1/32 8 0.1 Ymax*1.5]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    legend([PlotLAmp,PlotMAmp],'Location','southwest')
    set(PlotLAmp,...
        'DisplayName','L',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[204 0 0]/255)
    set(PlotMAmp,...
        'DisplayName','M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[119 172 48]/255)
%     set(PlotSAmp,...
%         'DisplayName','S',...
%         'LineWidth',2,...
%         'LineStyle','-',...
%         'Color',[0 0.45 0.74])

    subplot(4,3,10);
    PlotLPhase=semilogx(Freqs,abs(PlotL.Phase),'r');
    hold on;
    PlotMPhase=semilogx(Freqs,abs(PlotM.Phase),'g');
%     PlotSPhase=semilogx(Freqs,abs(PlotS.Phase),'b');
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Phase (deg)'); 
    axis([1/32 8 -20 200]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    %legend([PlotLPhase,PlotMPhase],'Location','northeast')
    set(PlotLPhase,...
        'DisplayName','L',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[204 0 0]/255)
    set(PlotMPhase,...
        'DisplayName','M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[119 172 48]/255)
%     set(PlotSPhase,...
%         'DisplayName','S',...
%         'LineWidth',2,...
%         'LineStyle','-',...
%         'Color',[0 0.45 0.74])

    subplot(4,3,[2 8]);
%     Ymax=max([PlotLMSum.Amplitude]);
%     Ymin=min([PlotLMSum.Amplitude]);
    if Ymax==0
        Ymax=10;
    end
    PlotTotAmp=loglog(Freqs,PlotLMSum.Amplitude,'r');
    hold on;
    %title({centerstrength;surroundstrength});
    %xlabel('Stimulus frequency (cpd)');
    ylabel('Response amplitude (a.u.)'); 
    axis([1/32 8 0.1 Ymax*1.5]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    legend(PlotTotAmp,'Location','southwest')
    set(PlotTotAmp,...
        'DisplayName','L+M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[.8 .8 .8])
    subplot(4,3,11);
    PlotTotPhase=semilogx(Freqs,abs(PlotLMSum.Phase),'r');
    hold on;
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Phase (deg)'); 
    axis([1/32 8 -20 200]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    %legend(PlotTotPhase,'Location','northeast')
    set(PlotTotPhase,...
        'DisplayName','L+M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[.8 .8 .8])

    subplot(4,3,[3 9]);
%     Ymax=max([PlotLMDiff.Amplitude]);
%     Ymin=min([PlotLMDiff.Amplitude]);
    if Ymax==0
        Ymax=10;
    end
    PlotLMDiffAmp=loglog(Freqs,PlotLMDiff.Amplitude,'r');
    hold on;
    %title({centerstrength;surroundstrength});
    %xlabel('Stimulus frequency (cpd)');
    ylabel('Response amplitude (a.u.)'); 
    axis([1/32 8 0.1 Ymax*1.5]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    legend(PlotLMDiffAmp,'Location','southwest')
    set(PlotLMDiffAmp,...
        'DisplayName','L-M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[.31 .31 .31])

    subplot(4,3,12);
    PlotLMDiffPhase=semilogx(Freqs,abs(PlotLMDiff.Phase),'r');
    hold on;
    %title({centerstrength;surroundstrength});
    xlabel('Stimulus frequency (cpd)');
    ylabel('Phase (deg)'); 
    axis([1/32 8 -20 200]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    set(gca,'XTick',[0.1 1 10]);
    %legend(PlotTotPhase,'Location','northeast')
    set(PlotLMDiffPhase,...
        'DisplayName','L-M',...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[.31 .31 .31])

    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,{centerstrength;surroundstrength;ori},'HorizontalAlignment','center','VerticalAlignment', 'top');
end
set(0,'DefaultFigureWindowStyle','default');