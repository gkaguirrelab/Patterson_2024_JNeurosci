function Cell=MatrixDoGfrom2D(Em)
%% IMPORT DATA
Cell=GenerateCell(Em);

cpAll=Cell.CenterWeights;
spAll=Cell.SurroundWeights;

%% MODEL UNIVERSE

Fs=512; %Sampling rate
numdeg=512;
x=-numdeg/2:1/Fs:numdeg/2-(1/Fs);
%Em=5; %Eccentricity in mm
Ed=(Em*1000)/200; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)
w=length(x);
Freq=(-w/2:(w/2)-1)*(Fs/w);

SigmasToRFRadius=3; %Let's say that ~3 SDs make up a RF radius
CSRadRatio=6; %Let's also say that a surround radius is ~6x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Finally, let's say that surround gain is ~3/4 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.3736*Em^2 + 0.01126*Em + 1.026);
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

DFRadius=0.002738*Em.^1.327; %DF radius data, in mm
DFRadiusDeg=((0.002738*Em.^1.327)*1000)/200; %This is in degrees.
CenterSigma=DFRadius/SigmasToRFRadius; %Here we assign a sigma for our cell given a particular eccentricity.
SurroundSigma=CSRadRatio*CenterSigma; %Ditto.

%ON or OFF cell? -1: OFF, +1: ON
CellType=1;

%% IMPORT DATA

AllCoords=Cell.Global.AllConeCoords(:,1:2);
LCoords=Cell.Global.LConeCoords;
MCoords=Cell.Global.MConeCoords;
LDoG=Cell.Global.DoG.LDoG;
MDoG=Cell.Global.DoG.MDoG;
TotalDoG=LDoG+MDoG;
AllCoordsDeg=(AllCoords*1000)/200;
LCoordsDeg=(LCoords*1000)/200;
MCoordsDeg=(MCoords*1000)/200;

%Round the mosaic locations so they can be assigned to a discrete matrix
RoundedCoords=ceil(AllCoordsDeg.*1000)/1000;
RoundedLCoords=ceil(LCoordsDeg.*1000)/1000;
RoundedMCoords=ceil(MCoordsDeg.*1000)/1000;

%% RF STRUCTURE

% Create matrix FXY into which to insert RoundedLocation values
Fx=zeros(10001);
lin=(-5:.001:5);
for f=1:length(Fx(:,1))
Fx(f,:)=lin;
end

Fy=zeros(10001);
lin=(-5:0.001:5);
for f=1:length(Fy(1,:))
Fy(:,f)=lin';
end
Fy=flipud(Fy(:,1:end));

FXY=zeros(10001,10001);

% Insert RoundedLocations values
for p=1:length(RoundedCoords)
    if isnan(RoundedCoords(p,1))==0
        testx(p,:)=find(round(Fx(1,:)-RoundedCoords(p,1),4)==0);
        testy(p,:)=find(round(Fy(:,1)-RoundedCoords(p,2),4)==0);
        FXY(testy(p,:),testx(p,:))=1;
    else
        testx(p,:)=0;
        testy(p,:)=0;
    end
end

%Add the DoG values associated with each pixel
DoGGauss=zeros(10001,10001);
for d=1:length(testx)
    if testx(d,:)>0
        DoGGauss(testy(d),testx(d))=TotalDoG(d);
    end
end
%Add the DoG values associated with each pixel
DoGGaussWithDiams=zeros(10001,10001);
for d=1:length(testx)
    if testx(d,:)>0
        DoGGaussWithDiams(testy(d)-6:testy(d)+6,testx(d)-6:testx(d)+6)=TotalDoG(d);
    end
end

%% PRESENT THE STIMULUS

Freqs=[1/2 1 2 4 8];
%SFLabels={'SF=1/128' 'SF=1/64' 'SF=1/32' 'SF=1/16' 'SF=1/8' 'SF=1/4' 'SF=1/2' 'SF=1' 'SF=2' 'SF=4' 'SF=8' 'SF=16' 'SF=32' 'SF=64' 'SF=128'};
SpatFreq=2;
theta = 135; %Grating orientation (degrees)
phase = 0; %Phase (degree)
phi = (phase/360)*2*pi;
thetaRad=(theta/360)*2*pi; % convert theta (orientation) to radians
X = linspace(-2.5,2.5,1001);
[Xm, Ym] = meshgrid(X, X);
Xt = Xm*cos(thetaRad); %Compute proportion of Xm for given orientation
Yt = Ym*sin(thetaRad); %Compute proportion of Ym for given orientation
XYt=[Xt+Yt]; %Sum X and Y components
grating = cos(2*pi*SpatFreq*XYt+phi);

for t=1:length(Freqs)
    SpatFreq=Freqs(t); %Spatial frequency (cpd) **THIS MUST BE A POWER OF 2**
    T=1/SpatFreq; %Fundamental period
    TotalPeriods=(max(X)-min(X))/T;
    A=1;
    L_Iso=A*cos(2*pi*SpatFreq*XYt+phi);
    %M_Iso=A*cos(2*pi*SpatFreq*XYt+phi);
    %S_Iso=A*cos(2*pi*SpatFreq*XYt+phi);
    
    %%%%%%%%%% OUTPUT RESPONSE %%%%%%%%%%%%%

    TotalLResp=conv2(DoGGauss,L_Iso,'same');
    %TotalMResp=conv2(MDoGGauss,M_Iso,'same');
    %TotalSResp=conv2(SDoGGauss,S_Iso,'same');
    %TotalResp=TotalLResp+TotalMResp;
    
    W=length(X);
    TotalXf=(fft(TotalLResp,W));
    TotalPowerXf=abs(TotalXf)/floor(W/2);
    TotalPhaseXf= angle(TotalXf)*180/pi;
    figure;plot(TotalPowerXf(1000:9000,5000));
end

W=length(X);
TotalXf=(fft(TotalLResp,W));
TotalPowerXf=abs(TotalXf)/floor(W/2);
TotalPhaseXf= angle(TotalXf)*180/pi;
[Amp, Index] = max(TotalPowerXf);
TotalData.Amplitude(t) = Amp;
TotalData.Frequency(t) = numdeg/2 + Freq(Index);
TotalData.Phase(t) = TotalPhaseXf(Index);