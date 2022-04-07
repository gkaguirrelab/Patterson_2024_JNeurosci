function Cell=Grid2DDoG_Selective(Cell,Universe,thetas,SpatFreqs,SFLines,OriCircles)
%% Import cell data

Em = Cell.Eccentricity;
Selectivity = 0.01;

LCenterCoords=Cell.Center.LCones.Coords;
LCenterCoordsDeg=(LCenterCoords*1000)/200;
LCenterCoordsDeg(isnan(LCenterCoordsDeg(:,1)),:) = [];
LCenterKernel=(~isnan(Cell.Center.LCones.Coords(:,1))).*(Cell.Center.CenterKernel);
LCenterKernel(LCenterKernel==0)=[];

MCenterCoords=Cell.Center.MCones.Coords;
MCenterCoordsDeg=(MCenterCoords*1000)/200;
MCenterCoordsDeg(isnan(MCenterCoordsDeg(:,1)),:) = [];
MCenterKernel=(~isnan(Cell.Center.MCones.Coords(:,1))).*(Cell.Center.CenterKernel);
MCenterKernel(MCenterKernel==0)=[];

LSurroundCoords=Cell.Surround.LCones.Coords;
LSurroundCoordsDeg=(LSurroundCoords*1000)/200;
LSurroundCoordsDeg(isnan(LSurroundCoordsDeg(:,1)),:) = [];
LSurroundKernel=(~isnan(Cell.Surround.LCones.Coords(:,1))).*(Cell.Surround.SurroundKernel);
LSurroundKernel(LSurroundKernel==0)=[];

MSurroundCoords=Cell.Surround.MCones.Coords;
MSurroundCoordsDeg=(MSurroundCoords*1000)/200;
MSurroundCoordsDeg(isnan(MSurroundCoordsDeg(:,1)),:) = [];
MSurroundKernel=(~isnan(Cell.Surround.MCones.Coords(:,1))).*(Cell.Surround.SurroundKernel);
MSurroundKernel(MSurroundKernel==0)=[];

CenterRadius=max(Cell.Center.CenterConeDistances);
CenterRadiusDeg=(CenterRadius*1000)/200;

SurroundRadius=max(Cell.Surround.SurroundConeDistances);
SurroundRadiusDeg=(SurroundRadius*1000)/200;

%% Shift selectivity
%Selectivity mechanism assesses the center weights compared to the surround
if Cell.Polarity > 0
    SelectivityGain = find(Cell.CenterWeights-Cell.SurroundWeights > 0); % 1 = L cones, 2 = M cones
    SelectivityLoss = find(Cell.CenterWeights-Cell.SurroundWeights < 0);
elseif Cell.Polarity < 0
    SelectivityGain = find(Cell.CenterWeights-Cell.SurroundWeights < 0);
    SelectivityLoss = find(Cell.CenterWeights-Cell.SurroundWeights > 0);
end

% %Selectivity mechanism just looks for the majority in the center alone
% if Cell.CenterWeights(1) > Cell.CenterWeights(2)
%     SelectivityGain = 1;
%     SelectivityLoss = 2;
% elseif Cell.CenterWeights(1) < Cell.CenterWeights(2)
%     SelectivityGain = 2;
%     SelectivityLoss = 1;
% end

Cell.SelectiveCenterWeights(SelectivityGain) = (1+Selectivity)*Cell.CenterWeights(SelectivityGain);
Cell.SelectiveCenterWeights(SelectivityLoss) = (1-Selectivity)*Cell.CenterWeights(SelectivityLoss);
Cell.SelectiveCenterWeights(3) = 0;

LCenterKernel = (Cell.SelectiveCenterWeights(1)-sum(LCenterKernel))/length(LCenterKernel)+LCenterKernel;
MCenterKernel = (Cell.SelectiveCenterWeights(2)-sum(MCenterKernel))/length(MCenterKernel)+MCenterKernel;

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
    sumAmp=LMSumpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    diffAmp=LMDiffpower(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    normAmp=(max(max(TotalLpower))+max(max(TotalMpower)))/100;

    lPhase=TotalLphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    mPhase=TotalMphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    sumPhase=LMSumphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);
    diffPhase=LMDiffphase(SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices);

    TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=lAmp/normAmp;
    TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=mAmp/normAmp;
    LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=sumAmp/normAmp;
    LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Amplitude=diffAmp/normAmp;

    TotalLData.(matlab.lang.makeValidName(OrientationLabel)).Phase=lPhase;%unwrap(lPhase);
    TotalMData.(matlab.lang.makeValidName(OrientationLabel)).Phase=mPhase;%unwrap(mPhase);
    LMSumData.(matlab.lang.makeValidName(OrientationLabel)).Phase=sumPhase;%unwrap(sumPhase);
    LMDiffData.(matlab.lang.makeValidName(OrientationLabel)).Phase=diffPhase;%unwrap(diffPhase);
end

for s=1:length(SpatFreqs)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    lAmp=TotalLpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mAmp=TotalMpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sumAmp=LMSumpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    diffAmp=LMDiffpower(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    normAmp=(max(max(TotalLpower))+max(max(TotalMpower)))/100;

    lPhase=TotalLphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    mPhase=TotalMphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    sumPhase=LMSumphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);
    diffPhase=LMDiffphase(OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices);

    TotalLData.(matlab.lang.makeValidName(SFLabel)).Amplitude=lAmp/normAmp;
    TotalMData.(matlab.lang.makeValidName(SFLabel)).Amplitude=mAmp/normAmp;
    LMSumData.(matlab.lang.makeValidName(SFLabel)).Amplitude=sumAmp/normAmp;
    LMDiffData.(matlab.lang.makeValidName(SFLabel)).Amplitude=diffAmp/normAmp;

    TotalLData.(matlab.lang.makeValidName(SFLabel)).Phase=lPhase;%unwrap(lPhase);
    TotalMData.(matlab.lang.makeValidName(SFLabel)).Phase=mPhase;%unwrap(mPhase);
    LMSumData.(matlab.lang.makeValidName(SFLabel)).Phase=sumPhase;%unwrap(sumPhase);
    LMDiffData.(matlab.lang.makeValidName(SFLabel)).Phase=diffPhase;%unwrap(diffPhase);
end

%% TEST FOR CHROMATIC OPPONENCY 

%Assess chromatic opponency using the phase relationship
AmpIndex=LMSumData.Theta0deg.Amplitude(1)-LMDiffData.Theta0deg.Amplitude(1);
PhaseIndex=abs(TotalLData.Theta0deg.Phase(1))-abs(TotalMData.Theta0deg.Phase(1));

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

% Cell.RFProfile.SResponse.CenterSigma=AllSc;
% Cell.RFProfile.SResponse.CenterAmplitude=cpS;
% Cell.RFProfile.SResponse.SurroundSigma=AllSs;
% Cell.RFProfile.SResponse.SurroundAmplitude=spS;

Cell.ResponseFunctions.LIsoResponse=TotalLData;
Cell.ResponseFunctions.MIsoResponse=TotalMData;
Cell.ResponseFunctions.LMSumResponse=LMSumData;
Cell.ResponseFunctions.LMDiffResponse=LMDiffData;
