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

%% CHOOSE DATASET SIZE AND ECCENTRICITY RANGE

FunctionDirectory=strcat('/Users/Lauren/Documents/MATLAB/Functions');
RunName=strcat('2DMidgetRun_ ',date);
DataDirectory=strcat('/Users/lauren/Desktop/OpponencyModel/Modeling Data/',RunName);
if (exist(DataDirectory) == 0)
    mkdir(DataDirectory);
end

NumberOfCells=1000; %Number of cells per group (5 groups total)
MinEccentricity=.5;
MaxEccentricity=10;

%% MAKE THE CELLS!
start=GetSecs;
for i=1:1
%     for f=1:NumberOfCells %NORMAL DISTRIBUTION
%         n=0;
%         while n==0;
%             g=normrnd(6.16,1.65,1,1); %Emulate the distribution of the empirical data
%             if g>0 && g<10
%                 n=1;
%                 Eccentricities(f,:)=g;
%             end
%         end
%     end
    Eccentricities=(MaxEccentricity-MinEccentricity).*rand(NumberOfCells,1) + MinEccentricity; %UNIFORM DISTRIBUTION
    CellNumber=linspace(1,length(Eccentricities),length(Eccentricities)); %How many cells per sample?

    for e=1:length(Eccentricities)
        EccentricityLabel=strcat('Cell',num2str(CellNumber(e)),'_Ecc',num2str(Eccentricities(e)),'mm');
        CellLabel=strcat('Cell',num2str(CellNumber(e)));
        Cell=Grid2DDoG_CrysS_Bipolar_H2(Eccentricities(e),Universe,thetas,SpatFreqs,SFLines,OriCircles);
        Cell.Eccentricity=Eccentricities(e);
        Data.(matlab.lang.makeValidName(EccentricityLabel))=Cell;
        disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
    end

    cd(DataDirectory);
    if (exist(strcat(RunName,'.mat')))==2
        %Data=load(strcat(RunName,'.mat')); %COMMENT OUT  FOR SPEED
        save(RunName,'-struct','Data','-append'); %COMMENT OUT  FOR SPEED
    elseif (exist(strcat(RunName,'.mat')))==0
        save(RunName,'-struct','Data');
    end

    clearvars Data Eccentricities;
    disp(strcat('Group ',num2str(i),' complete.'));
    cd(FunctionDirectory);
end
stop=GetSecs;
disp(strcat('Finished in .',num2str((stop-start)/3600),' hours. Good job!'));
%% 
% figure;
% xlabel('Center L/(L+M)');
% ylabel('Eccentricity'); 
% zlabel('Surround L/(L+M)'); 
% axis([0 1 0.25 10 0 1]);
% view([45 22]);
% grid on;
% hold on;
% Unity=fill3([0,1,1,0],[0,0,10,10],[0,1,1,0],'k');
% set(Unity,...
%     'LineStyle','none',...
%     'FaceAlpha',0.125)
% for e=1:length(Eccentricities)
%     EccentricityLabel=strcat('Ecc',num2str(Eccentricities(e)),'mm');
%     CellLabel=strcat('Cell',num2str(CellNumber(e)));
%     PlotEccentricity=Eccentricities(e);
%     PlotCenterWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).CenterWeights(1);
%     PlotSurroundWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).SurroundWeights(1);
%     PlotChromTag=Data.(matlab.lang.makeValidName(EccentricityLabel)).ChromTag;
%     if strcmp(PlotChromTag,'L-dominated')==1
%         PlotL=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'ro');
%         set(PlotL,...
%             'LineStyle','none',...
%             'Color','k',...
%             'MarkerFaceColor',[0.64 0.08 0.18],...
%             'MarkerSize',8,...
%             'Marker','o')
%     elseif strcmp(PlotChromTag,'M-dominated')==1
%         PlotM=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'go');
%         set(PlotM,...
%             'LineStyle','none',...
%             'Color','k',...
%             'MarkerFaceColor',[0.47 0.67 0.19],...
%             'MarkerSize',8,...
%             'Marker','o')
%     elseif strcmp(PlotChromTag,'Achromatic')==1
%         PlotA=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'ko');
%         set(PlotA,...
%             'LineStyle','none',...
%             'Color','k',...
%             'MarkerFaceColor',[1 0.84 0],...
%             'MarkerSize',8,...
%             'Marker','o')
%     end
% end
% 
