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

%% GENERATE GRATING

%Freqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8];
Freqs=[1];
Thetas=[0];
Phases=0:10:350;

for f=1:length(Freqs)
    SpatFreq=Freqs(f);
    SFLabel=strcat('SF',num2str(SpatFreq),'cpd');
    for t=1:length(Thetas)
        theta = Thetas(t); %Grating orientation (degrees)
        OriLabel=strcat('Ori',num2str(theta),'deg');
        for p=1:length(Phases)
            phase = Phases(p); %Phase (degree)
            PhaseLabel=strcat('Phase',num2str(phase),'deg');
            phi = (phase/360)*2*pi;
            thetaRad=(theta/360)*2*pi; % convert theta (orientation) to radians
            [Xm, Ym] = meshgrid(X, X);
            Xt = Xm*cos(thetaRad); %Compute proportion of Xm for given orientation
            Yt = Ym*sin(thetaRad); %Compute proportion of Ym for given orientation
            XYt=Xt+Yt; %Sum X and Y components
            grating=cos(2*pi*SpatFreq*XYt+phi);
            Stimulus.(matlab.lang.makeValidName(SFLabel)).(matlab.lang.makeValidName(OriLabel)).Grating(:,:,p)=grating;
        end
    end
end

%% COMPUTE FFT FOR EACH GRATING (SF x ORI x THETA x PHASE)
for f=1:length(Freqs)
    SpatFreq=Freqs(f);
    SFLabel=strcat('SF',num2str(SpatFreq),'cpd');
    for t=1:length(Thetas)
        theta = Thetas(t); %Grating orientation (degrees)
        OriLabel=strcat('Ori',num2str(theta),'deg');
        for p=1:length(Phases)
            grating=Stimulus.(matlab.lang.makeValidName(SFLabel)).(matlab.lang.makeValidName(OriLabel)).Grating(:,:,p);
            gFFT=fftshift(fft2(ifftshift(grating)));
            Stimulus.(matlab.lang.makeValidName(SFLabel)).(matlab.lang.makeValidName(OriLabel)).gFFT(:,:,p)=gFFT;
        end
    end
end

%% PAIRWISE PRODUCT WITH THE RELEVANT RF FFTS (this should eventually go in the actual DoG script)

LgFFTprod=gFFT.*TotalLFFT;
LgPower=abs(LgFFTprod);
LgPhase=angle(LgFFTprod)*180/pi;

LgMaxInd=find(LgPower==(max(max(LgPower))));
[LgMaxI1, LgMaxJ1]=ind2sub(size(LgPower),LgMaxInd(1));
[LgMaxI2, LgMaxJ2]=ind2sub(size(LgPower),LgMaxInd(2));

MgFFTprod=gFFT.*TotalMFFT;
MgPower=abs(MgFFTprod);
MgPhase=angle(MgFFTprod)*180/pi;

MgMaxInd=find(MgPower==(max(max(MgPower))));
[MgMaxI1, MgMaxJ1]=ind2sub(size(MgPower),MgMaxInd(1));
[MgMaxI2, MgMaxJ2]=ind2sub(size(MgPower),MgMaxInd(2));
