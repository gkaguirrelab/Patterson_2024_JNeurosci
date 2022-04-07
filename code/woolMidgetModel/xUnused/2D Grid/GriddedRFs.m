%function GriddedRFs(Em)
%% MODEL UNIVERSE

NumDegs=5; %Number of degrees across a dimension (to create an NumDeg-by-NumDeg patch of retina)
XRes=1000; %How resolved is the patch? (i.e., matrix length = XRes)
Fs=XRes/NumDegs; %Sampling frequency
X = -NumDegs/2:1/Fs:NumDegs/2-(1/Fs); %Where in space is this NumDeg-by-NumDeg patch?
W=length(X); %Determine length of X dimension for FFT frequency determination
Freq=(-W/2:(W/2)-1)*(Fs/W); %Frequency domain values for FFT

%% IMPORT DATA
%Cell=GenerateBipolarLayer(Em);

AllCoords=Cell.Global.AllConeCoords(:,1:2);
LCoords=Cell.Global.LConeCoords;
MCoords=Cell.Global.MConeCoords;
LCenterCoords=Cell.Center.LCones.Coords;
MCenterCoords=Cell.Center.MCones.Coords;
LSurroundCoords=Cell.Surround.LCones.Coords;
MSurroundCoords=Cell.Surround.MCones.Coords;
LCenterKernel=(~isnan(Cell.Center.LCones.Coords(:,1))).*(Cell.Center.CenterKernel);
MCenterKernel=(~isnan(Cell.Center.MCones.Coords(:,1))).*(Cell.Center.CenterKernel);
LSurroundKernel=(~isnan(Cell.Surround.LCones.Coords(:,1))).*(Cell.Surround.SurroundKernel);
MSurroundKernel=(~isnan(Cell.Surround.MCones.Coords(:,1))).*(Cell.Surround.SurroundKernel);
LDoG=Cell.Global.DoG.LDoG;
MDoG=Cell.Global.DoG.MDoG;
TotalDoG=LDoG+MDoG;
AllCoordsDeg=(AllCoords*1000)/200;
LCoordsDeg=(LCoords*1000)/200;
MCoordsDeg=(MCoords*1000)/200;
LCenterCoordsDeg=(LCenterCoords*1000)/200;
MCenterCoordsDeg=(MCenterCoords*1000)/200;
LSurroundCoordsDeg=(LSurroundCoords*1000)/200;
MSurroundCoordsDeg=(MSurroundCoords*1000)/200;

%% CLEAN UP VECTORS
%Clean out the NaNs from the coordinate arrays
LCoordsDeg(isnan(LCoordsDeg(:,1)),:) = [];
MCoordsDeg(isnan(MCoordsDeg(:,1)),:) = [];

%Fix up the kernels and get rid of extra zeros
[~,v]=max(cumsum(LDoG'~=0));
for h=1:v
    LDoG(LDoG==0)=[];
end
LDoG=padarray(LDoG,length(LCoordsDeg)-length(LDoG),'post');
[~,v]=max(cumsum(MDoG'~=0));
for h=1:v
    MDoG(MDoG==0)=[];
end
MDoG=padarray(MDoG,length(MCoordsDeg)-length(MDoG),'post');

[~,v]=max(cumsum(LCenterKernel'~=0));
for h=1:v
    LCenterKernel(LCenterKernel==0)=[];
end
LCenterKernel=padarray(LCenterKernel,length(LCoordsDeg)-length(LCenterKernel),'post');
[~,v]=max(cumsum(MCenterKernel'~=0));
for h=1:v
    MCenterKernel(MCenterKernel==0)=[];
end
MCenterKernel=padarray(MCenterKernel,length(MCoordsDeg)-length(MCenterKernel),'post');

[~,v]=max(cumsum(LSurroundKernel'~=0));
for h=1:v
    LSurroundKernel(LSurroundKernel==0)=[];
end
LSurroundKernel=padarray(LSurroundKernel,length(LCoordsDeg)-length(LSurroundKernel),'post');
[~,v]=max(cumsum(MSurroundKernel'~=0));
for h=1:v
    MSurroundKernel(MSurroundKernel==0)=[];
end
MSurroundKernel=padarray(MSurroundKernel,length(MCoordsDeg)-length(MSurroundKernel),'post');

%% GENERATE L AND M RESPONSE SURFACES
%Assign a meshgrid to interpolate the data over
[xq,yq] = meshgrid(X,X);

%Center only
GriddedLCenterInputs=griddata(LCoordsDeg(:,1),LCoordsDeg(:,2),LCenterKernel,xq,yq);
GriddedLCenterInputs(isnan(GriddedLCenterInputs))=0;
GriddedMCenterInputs=griddata(MCoordsDeg(:,1),MCoordsDeg(:,2),MCenterKernel,xq,yq);
GriddedMCenterInputs(isnan(GriddedMCenterInputs))=0;

%Surround only
GriddedLSurroundInputs=griddata(LCoordsDeg(:,1),LCoordsDeg(:,2),LSurroundKernel,xq,yq);
GriddedLSurroundInputs(isnan(GriddedLSurroundInputs))=0;
GriddedMSurroundInputs=griddata(MCoordsDeg(:,1),MCoordsDeg(:,2),MSurroundKernel,xq,yq);
GriddedMSurroundInputs(isnan(GriddedMSurroundInputs))=0;

%% MAKE A STIMULUS GRATING

Freqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32 64 128];
%SFLabel=strcat('SF=',num2str(Freqs(t)));
t=8;
%for t=1:length(Freqs)
SpatFreq=Freqs(t);
SpatFreq1=1;
SpatFreq2=2;
theta = 30; %Grating orientation (degrees)
phase = 90; %Phase (degree)
phi = (phase/360)*2*pi;
thetaRad=(theta/360)*2*pi; % convert theta (orientation) to radians
[Xm, Ym] = meshgrid(X, X);
Xt = Xm*cos(thetaRad); %Compute proportion of Xm for given orientation
Yt = Ym*sin(thetaRad); %Compute proportion of Ym for given orientation
XYt=[Xt+Yt]; %Sum X and Y components
grating = cos(2*pi*SpatFreq*XYt+phi);
% grating1=cos(2*pi*SpatFreq1*XYt+phi);
% grating2=cos(2*pi*SpatFreq2*XYt+phi);
% grating=grating1-grating2;

%Plot it if you want
figure;imagesc(X,X,grating);axis square;

%% CONVOLVE CENTERS AND SURROUNDS

% LcResp=ifft2(fft2(GriddedLCenterInputs).*fft2(grating));
% McResp=ifft2(fft2(GriddedMCenterInputs).*fft2(grating));
% LsResp=ifft2(fft2(GriddedLSurroundInputs).*fft2(grating));
% MsResp=ifft2(fft2(GriddedMSurroundInputs).*fft2(grating));

LcResp=conv_fft2(GriddedLCenterInputs,grating,'same');
McResp=conv_fft2(GriddedMCenterInputs,grating,'same');
LsResp=conv_fft2(GriddedLSurroundInputs,grating,'same');
MsResp=conv_fft2(GriddedMSurroundInputs,grating,'same');

TotalLResp=LcResp-.75*LsResp;
TotalMResp=McResp-.75*MsResp;
%TotalSResp=ScResp+SsResp;
%TotalResp=TotalLResp+TotalMResp+TotalSResp;
%TotalResp=TotalLResp+TotalMResp;
LMDiff=TotalLResp-TotalMResp;

%% FFT

TotalLFFT=fftshift(fft2(TotalLResp));
TotalLPower=abs(TotalLFFT)/floor(W/2);
TotalLPhase=unwrap(angle(TotalLFFT))*180/pi;
[Yi,Xi]=find(TotalLPower==(max(max(TotalLPower))));
Index=[Yi(1),Xi(1)];
Amp=TotalLPower(Index(1),Index(2));
Frequency=sqrt(Freq(Index(1))^2+Freq(Index(2))^2);
Angle=atan(Freq(Index(1))/Freq(Index(2)))*(180/pi);
if Angle<0
    Angle=180+Angle;
end
TotalLData.Amplitude(t) = Amp;
TotalLData.Frequency(t) = Frequency;
TotalLData.Phase(t) = TotalLPhase(Index(1),Index(2));

TotalMFFT=fftshift(fft2(TotalMResp));
TotalMPower=abs(TotalMFFT);
TotalMPhase=angle(TotalMFFT)*180/pi;
[Yi,Xi]=find(TotalMPower==(max(max(TotalMPower))));
Index=[Yi(1),Xi(1)];
Amp=TotalMPower(Index(1),Index(2));
Frequency=sqrt(Freq(Index(1))^2+Freq(Index(2))^2);
Angle=atan(Freq(Index(1))/Freq(Index(2)))*(180/pi);
if Angle<0
    Angle=180+Angle;
end
TotalMData.Amplitude(t) = Amp;
TotalMData.Frequency(t) = Frequency;
TotalMData.Phase(t) = TotalMPhase(Index(1),Index(2));
%end

%% IDEALIZED RF
CSigma=max([LCenterCoordsDeg, MCenterCoordsDeg]);
CGauss=exp(-(((Xm.^2)+(Ym.^2))./(2*CSigma^2))); % formula for 2D gaussian
CGauss=CGauss/(sum(sum(CGauss)));
SSigma=max([LSurroundCoordsDeg, MSurroundCoordsDeg]);
SGauss=exp(-(((Xm.^2)+(Ym.^2))./(2*SSigma^2))); % formula for 2D gaussian
SGauss=(SGauss/(sum(sum(SGauss))))*(-0.75);
CResp=conv_fft2(CGauss,grating,'same');
SResp=conv_fft2(SGauss,grating,'same');
TotalResp=CResp;
TotalFFT=fftshift(fft2(TotalResp));
TotalPower=abs(TotalFFT);
TotalPhase=angle(TotalFFT)*180/pi;
[Yi,Xi]=find(TotalPower==(max(max(TotalPower))));
Index=[Yi(1),Xi(1)];
Amp=TotalPower(Index(1),Index(2));
Frequency=sqrt(Freq(Index(1))^2+Freq(Index(2))^2);
Angle=atan(Freq(Index(1))/Freq(Index(2)))*(180/pi);
if Angle<0
    Angle=180+Angle;
end
TotalData.Amplitude(t) = Amp;
TotalData.Frequency(t) = Frequency;
TotalData.Phase(t) = TotalPhase(Index(1),Index(2));

%% PLOT THE RESULT
figure;
subplot(2,3,1);
imagesc(X,X,GriddedLInputs)
axis square
axis([-2.5 2.5 -2.5 2.5]);
xlabel('x (deg)');
ylabel('y (deg)');
title('L RESPONSE PROFILE');
subplot(2,3,4);
imagesc(X,X,GriddedMInputs)
axis square
axis([-2.5 2.5 -2.5 2.5]);
xlabel('x (deg)');
ylabel('y (deg)');
title('M RESPONSE PROFILE');

subplot(2,3,2);
imagesc(X,X,LWithGrating)
axis square
axis([-2.5 2.5 -2.5 2.5]);
xlabel('x (deg)');
ylabel('y (deg)');
title('L RESPONSE');
subplot(2,3,5);
imagesc(X,X,MWithGrating)
axis square
axis([-2.5 2.5 -2.5 2.5]);
xlabel('x (deg)');
ylabel('y (deg)');
title('M RESPONSE');

subplot(2,3,3);
imagesc(X,X,sLFFTpower)
axis square
title('L FFT');
subplot(2,3,6);
imagesc(X,X,sMFFTpower)
axis square
title('M FFT');

figure;
subplot(1,2,1)
surf(X,X,LWithGrating);
shading interp;
axis square;
axis([-2.5 2.5 -2.5 2.5 -.05 0.125]);
xlabel('x (deg)');
ylabel('y (deg)');
zlabel('Normalized amplitude');
title('L RESPONSE');
subplot(1,2,2)
surf(X,X,MWithGrating);
shading interp;
axis square;
axis([-2.5 2.5 -2.5 2.5 -.05 0.125]);
xlabel('x (deg)');
ylabel('y (deg)');
zlabel('Normalized amplitude');
title('M RESPONSE');

%% PLOT SF TUNING
PlotCell=Output;
CellLAmpData=PlotCell.L.Amp;
CellMAmpData=PlotCell.M.Amp;
CellLPhaseData=PlotCell.L.Phase;
CellMPhaseData=PlotCell.M.Phase;

Freqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32 64 128];
SFLabels={'SF=1/128' 'SF=1/64' 'SF=1/32' 'SF=1/16' 'SF=1/8' 'SF=1/4' 'SF=1/2' 'SF=1' 'SF=2' 'SF=4' 'SF=8' 'SF=16' 'SF=32' 'SF=64' 'SF=128'};

figure;
subplot(2,1,1);
Ymax=max([CellLAmpData CellMAmpData]);
if Ymax==0
    Ymax=10;
end
PlotLAmp=loglog(Freqs,CellLAmpData,'ro');
hold on;
PlotMAmp=loglog(Freqs,CellMAmpData,'go');
%PlotSAmp=loglog(Freqs,TotalSData.Amplitude,'bo');
%title({centerstrength;surroundstrength});
xlabel('Stimulus frequency (cpd)');
ylabel('Response amplitude (a.u.)'); 
axis([0.005 150 50 200]);
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
% set(PlotSAmp,...
%     'DisplayName','S',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[0 0.45 0.74])

subplot(2,1,2);
PlotLPhase=semilogx(Freqs,abs(CellLPhaseData),'ro');
hold on;
PlotMPhase=semilogx(Freqs,abs(CellMPhaseData),'go');
% PlotSPhase=semilogx(Freqs,abs(TotalSData.Phase),'go');
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
% set(PlotSPhase,...
%     'DisplayName','S',...
%     'LineWidth',.5,...
%     'LineStyle',':',...
%     'Color','k',...
%     'MarkerFaceColor',[0 0.45 0.74])

%% WORKSHOP
% 
% LFFT=fft2(LWithGrating);
% sLFFT=fftshift(LFFT);
% sLFFTpower=abs(sLFFT);
% [lAmp, lIndex]=max(max(sLFFTpower));
% sLFFTphase=angle(sLFFT)*180/pi;
% Output.L.Amp=lAmp;
% Output.L.Phase=sLFFTphase(lIndex);
% 
% MFFT=fft2(MWithGrating);
% sMFFT=fftshift(MFFT);
% sMFFTpower=abs(sMFFT);
% [mAmp, mIndex]=max(max(sMFFTpower));
% sMFFTphase=angle(sMFFT)*180/pi;
% Output.M.Amp=mAmp;
% Output.M.Phase=sMFFTphase(mIndex);
% 
% StimFFT=fft2(grating);
% sStimFFT=fftshift(StimFFT);
% 
% LcInputsFFT=fft2(GriddedLCenterInputs);
% sLcInputsFFT=fftshift(LcInputsFFT);
% McInputsFFT=fft2(GriddedMCenterInputs);
% sMcInputsFFT=fftshift(McInputsFFT);
% 
% LsInputsFFT=fft2(GriddedLSurroundInputs);
% sLsInputsFFT=fftshift(LsInputsFFT);
% MsInputsFFT=fft2(GriddedMSurroundInputs);
% sMsInputsFFT=fftshift(MsInputsFFT);
% 
% LcResp=sLcInputsFFT.*sStimFFT;
% sLcRespPower=abs(LcResp);
% [lcAmp, lcIndex]=max(max(sLcRespPower));
% sLcRespPhase=angle(LcResp)*180/pi;
% Output.L.Center.Amp=lcAmp;
% Output.L.Center.Phase=sLcRespPhase(lcIndex);
% McResp=sMcInputsFFT.*sStimFFT;
% sMcRespPower=abs(McResp);
% [mcAmp, mcIndex]=max(max(sMcRespPower));
% sMcRespPhase=angle(McResp)*180/pi;
% Output.M.Center.Amp=mcAmp;
% Output.M.Center.Phase=sMcRespPhase(mcIndex);
% 
% LsResp=sLsInputsFFT.*sStimFFT;
% sLsRespPower=abs(LsResp);
% [lcAmp, lcIndex]=max(max(sLsRespPower));
% sLsRespPhase=angle(LsResp)*180/pi;
% Output.L.Surround.Amp=lcAmp;
% Output.L.Surround.Phase=sLsRespPhase(lcIndex);
% MsResp=sMsInputsFFT.*sStimFFT;
% sMsRespPower=abs(MsResp);
% [mcAmp, mcIndex]=max(max(sMsRespPower));
% sMsRespPhase=angle(MsResp)*180/pi;
% Output.M.Surround.Amp=mcAmp;
% Output.M.Surround.Phase=sMsRespPhase(mcIndex);

% LcResp=1/Fs*conv2(GriddedLCenterInputs,grating,'same');
% McResp=1/Fs*conv2(GriddedMCenterInputs,grating,'same');
% %ScResp=1/Fs*conv2(Sc,S_Iso,'same');
% 
% LsResp=1/Fs*conv2(GriddedLSurroundInputs,grating,'same');
% MsResp=1/Fs*conv2(GriddedMSurroundInputs,grating,'same');
% %SsResp=1/Fs*conv2(Ss,S_Iso,'same');
% 
% %% MULTIPLY SURFACES WITH THE GRATING
% 
% %Total
% LWithGrating=grating.*GriddedLInputs;
% LWithGrating(isnan(LWithGrating))=0;
% MWithGrating=grating.*GriddedMInputs;
% MWithGrating(isnan(MWithGrating))=0;
% 
% %Center only
% LCenterWithGrating=grating.*GriddedLCenterInputs;
% LCenterWithGrating(isnan(LCenterWithGrating))=0;
% MCenterWithGrating=grating.*GriddedMCenterInputs;
% MCenterWithGrating(isnan(MCenterWithGrating))=0;
% 
% %Surround only
% LSurroundWithGrating=grating.*GriddedLSurroundInputs;
% LSurroundWithGrating(isnan(LSurroundWithGrating))=0;
% MSurroundWithGrating=grating.*GriddedMSurroundInputs;
% MSurroundWithGrating(isnan(MSurroundWithGrating))=0;

% %Total
% GriddedLInputs=griddata(LCoordsDeg(:,1),LCoordsDeg(:,2),LDoG,xq,yq);
% GriddedLInputs(isnan(GriddedLInputs))=0;
% GriddedMInputs=griddata(MCoordsDeg(:,1),MCoordsDeg(:,2),MDoG,xq,yq);
% GriddedMInputs(isnan(GriddedMInputs))=0;
