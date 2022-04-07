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
%thetas=[0 30 60 90 120 150];
SpatFreqs=[1/128 1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8];

%thetas=[0 45 90 135];
thetas=[0];
SFLines=ComputeSFLine(thetas,W,Freq);
OriCircles=ComputeOriCircle(SpatFreqs,W,Freq);
%%
Data = load('/Users/lauren/Documents/SCHOOL/Dacey/Midget Paper/FieldComparison/FieldComparison.mat');
Cells=fieldnames(Data);

for i=1:2%length(Cells)
    Cell=Data.(matlab.lang.makeValidName(Cells{i}));
    EccentricityLabel=strcat(Cells{i});
    CellLabel=strcat('Cell',num2str(i));
    Cell=Grid2DDoG_Selective(Cell,Universe,thetas,SpatFreqs,SFLines,OriCircles);

    SelectiveData.(matlab.lang.makeValidName(EccentricityLabel))=Cell;
    disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
end