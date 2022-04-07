function Parasol=GenerateParasol_BipolarLayer(Em)
%% INITIATE PARASOL CELL SPECS
%Identify the size of the parasol cell
DFRadius=-0.000565*(Em^2)+(0.01863*Em)+0.00146;

%% GENERATE THE MOSAICS
ConeMosaic=GenerateMosaic(Em);
BipolarMosaic=GenerateBipolarMosaic(Em);

%% BIPOLAR INPUTS TO GANGLION CELL

%Choose a central subset of the cone mosaic locations as your 
%substrate for the bipolar locations. This keeps processing speeds short...
BipolarSubstrate=BipolarMosaic.AllCoords(1:500,:);

%Determine the number of bipolars per ganglion cell
BipolarToGanglionIndex=find(BipolarSubstrate(:,3)<DFRadius);
BipolarCells=BipolarSubstrate(find(BipolarSubstrate(:,3)<DFRadius),1:2);
NumberofBipolarInputs=length(BipolarCells);

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Location=BipolarCells(b,:);
end

%% CONE INPUTS TO EACH BIPOLAR

%General specs
CSRadRatio=15;
CSStrengthRatio=.75;
ConesToBipolarCenter=7;
ConesToBipolarSurround=(CSRadRatio^2)*ConesToBipolarCenter;

%Compute the distances of each cone from each bipolar and choose the ones
%that make up the inputs to that particular bipolar
MosaicL=~isnan(ConeMosaic.LConeCoords(:,1));
MosaicM=~isnan(ConeMosaic.MConeCoords(:,1));
ConeSubstrate=ConeMosaic.AllConeCoords(1:5000,:);
LCones=MosaicL(1:length(ConeSubstrate));
MCones=MosaicM(1:length(ConeSubstrate));

for d=1:length(BipolarCells)
    BipolarDistances(:,d)=pdist2(ConeSubstrate(:,1:2),[BipolarCells(d,1),BipolarCells(d,2)]);
    [SortDistances(:,d), SortIndices(:,d)]=sort(BipolarDistances(:,d));
    BipolarConeCenterIndices(:,d)=SortIndices(1:ConesToBipolarCenter,d);
    BipolarConeSurroundIndices(:,d)=SortIndices(1:ConesToBipolarSurround,d);
end

%Determine the L and M cones for each bipolar
for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Total.Locations=ConeSubstrate(BipolarConeCenterIndices(:,b),:);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,1)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,2)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,1)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,2)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Total.Locations=ConeSubstrate(BipolarConeSurroundIndices(:,b),:);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,2)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,2)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
end

%Weight the cone inputs to each bipolar
for x=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(x));
    CenterX=SortDistances(1:ConesToBipolarCenter,x);
    SurroundX=SortDistances(1:ConesToBipolarSurround,x);
    
    %Assign a sigma to determine the shape of the kernels
    cSigma=max(CenterX)/1; %(we assume ~1 SDs make an RF radius)
    sSigma=max(SurroundX)/1;

    %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
    %from the origin...
    CenterKernel(:,x)=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
    NormalizedCenterKernel(:,x)=CenterKernel/(sum(CenterKernel));
    SurroundKernel(:,x)=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
    NormalizedSurroundKernel(:,x)=(SurroundKernel/(sum(SurroundKernel)));
    
    %Save the kernels to substruct
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Kernel=NormalizedCenterKernel(:,x);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Kernel=NormalizedSurroundKernel(:,x);
end

%Compute the DoG for each bipolar
PaddedNormalizedCenterKernel=padarray(NormalizedCenterKernel,size(NormalizedSurroundKernel,1)-size(NormalizedCenterKernel,1),'post');
DoG=PaddedNormalizedCenterKernel-(CSStrengthRatio*NormalizedSurroundKernel);
    
%Assign cones to the kernel (based on location) to determine the weight of each subtype
for h=1:NumberofBipolarInputs
    CenterLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*DoG(:,h);
    CenterMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*DoG(:,h);  
end

for z=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(z));
    %Save the DoG to substruct
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).DoG.Kernel=DoG(:,z);
    
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.LWeight=CenterLWeight(:,z);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.LWeight=SurroundLWeight(:,z);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).DoG.LCones.LWeight=TotalLWeight(:,z);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.MWeight=CenterMWeight(:,z);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.MWeight=SurroundMWeight(:,z);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).DoG.MCones.MWeight=TotalMWeight(:,z);
end

%Compute the discrete number of cells in the center and the surround (since
%there's overlap in each bipolar)
UniqueCenterIndices=unique(BipolarConeCenterIndices);
UniqueSurroundIndices=unique(BipolarConeSurroundIndices);
UniqueLCenterIndices=LCones(1:length(UniqueCenterIndices)).*UniqueCenterIndices;
UniqueLCenterIndices(UniqueLCenterIndices==0)=[];
UniqueMCenterIndices=MCones(1:length(UniqueCenterIndices)).*UniqueCenterIndices;
UniqueMCenterIndices(UniqueMCenterIndices==0)=[];
UniqueLSurroundIndices=LCones(1:length(UniqueSurroundIndices)).*UniqueSurroundIndices;
UniqueLSurroundIndices(UniqueLSurroundIndices==0)=[];
UniqueMSurroundIndices=MCones(1:length(UniqueSurroundIndices)).*UniqueSurroundIndices;
UniqueMSurroundIndices(UniqueMSurroundIndices==0)=[];

UniqueCenterCones=ConeSubstrate(UniqueCenterIndices,:);
UniqueSurroundCones=ConeSubstrate(UniqueSurroundIndices,:);
UniqueLCenterCones=ConeSubstrate(UniqueLCenterIndices,:);
UniqueMCenterCones=ConeSubstrate(UniqueMCenterIndices,:);
UniqueLSurroundCones=ConeSubstrate(UniqueLSurroundIndices,:);
UniqueMSurroundCones=ConeSubstrate(UniqueMSurroundIndices,:);

cnL=length(UniqueLCenterIndices);
cnM=length(UniqueMCenterIndices);
cnS=0;
snL=length(UniqueLSurroundIndices);
snM=length(UniqueMSurroundIndices);
snS=0;

%% WEIGHT BIPOLAR INPUTS TO THE PARASOL

%Apply Gaussian weight to each bipolar cell
BCdist=pdist2(BipolarCells,[0,0]);
BCsigma=max(BCdist)/1;
BCKernel=exp((-1/(2*BCsigma^2))*(abs(BCdist).^2));

%Apply the weight of each bipolar to the respective weights of each cone
%input
for g=1:length(BipolarCells)
    ADJCenterL(:,g)=CenterLWeight(:,g)*BCKernel(g);
    ADJSurroundL(:,g)=SurroundLWeight(:,g)*BCKernel(g);
    ADJCenterM(:,g)=CenterMWeight(:,g)*BCKernel(g);
    ADJSurroundM(:,g)=SurroundMWeight(:,g)*BCKernel(g);
end

TotalLCenter=sum(sum(ADJCenterL));
TotalLSurround=sum(sum(ADJSurroundL));
TotalMCenter=sum(sum(ADJCenterM));
TotalMSurround=sum(sum(ADJSurroundM));
TotalCenter=TotalLCenter+TotalMCenter;
TotalSurround=TotalLSurround+TotalMSurround;

cpL=TotalLCenter/TotalCenter;
cpM=TotalMCenter/TotalCenter;
cpS=0;
spL=TotalLSurround/TotalSurround;
spM=TotalMSurround/TotalSurround;
spS=0;

cnAll=[cnL cnM cnS];
snAll=[snL snM snS];
cpAll=[cpL cpM cpS];
spAll=[spL spM spS];

% Use these values for determining the center and surround sigmas for the
% DoG. They need to be converted to deg
MaxCenterDist=max(UniqueCenterCones(:,3));
MaxSurroundDist=max(UniqueSurroundCones(:,3));
MaxCenterDistDeg=(MaxCenterDist*1000)/200;
MaxSurroundDistDeg=(MaxSurroundDist*1000)/200;

%% SAVE TO OUTPUT

% GLOBAL

Parasol.ConeLayer.Mosaic=ConeMosaic;
Parasol.BipolarLayer=BipolarLayer;
Parasol.BipolarLayer.BipolarWeights=BCKernel;
Parasol.BipolarLayer.Mosaic=BipolarMosaic;

Parasol.CenterConeCount=[cnL cnM cnS];
Parasol.SurroundConeCount=[snL snM snS];
Parasol.CenterWeights=cpAll;
Parasol.SurroundWeights=spAll;

% CENTER
Parasol.Center.Radius=MaxCenterDist;
Parasol.Center.ConesToCenter=length(UniqueCenterIndices);
Parasol.Center.CenterConeCoords=UniqueCenterCones;
Parasol.Center.LCones.Count=cnL;
Parasol.Center.LCones.Coords=UniqueLCenterCones;
Parasol.Center.LCones.TotalWeight=cpL;
Parasol.Center.MCones.Count=cnM;
Parasol.Center.MCones.Coords=UniqueMCenterCones;
Parasol.Center.MCones.TotalWeight=cpM;
Parasol.Center.SCones.Count=cnS;
Parasol.Center.SCones.TotalWeight=cpS;

% SURROUND
Parasol.Surround.Radius=MaxSurroundDist;
Parasol.Surround.ConesToSurround=length(UniqueSurroundIndices);
Parasol.Surround.SurroundConeCoords=UniqueSurroundCones;
Parasol.Surround.LCones.Count=cnL;
Parasol.Surround.LCones.Coords=UniqueLSurroundCones;
Parasol.Surround.LCones.TotalWeight=cpL;
Parasol.Surround.MCones.Count=cnM;
Parasol.Surround.MCones.Coords=UniqueMSurroundCones;
Parasol.Surround.MCones.TotalWeight=cpM;
Parasol.Surround.SCones.Count=cnS;
Parasol.Surround.SCones.TotalWeight=cpS;

%% PLOT THE DATA
% 
% %Plot the 1-mm2 patch and its L,M,S cones
% figure;
% hold on;
% MapL=plot(LCoords(:,1),LCoords(:,2),'r.');
% MapM=plot(MCoords(:,1),MCoords(:,2),'g.');
% MapS=plot(SCoords(:,1),SCoords(:,2),'b.');
% axis square
% set(gca, 'color', [0.25 0.25 0.25]);
% set(MapL,...
%     'DisplayName','L Cones',...
%     'MarkerFaceColor',[.64 .08 .18],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',3,...
%     'Marker','o')
% set(MapM,...
%     'DisplayName','M Cones',...
%     'MarkerFaceColor',[.47 .67 .19],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',3,...
%     'Marker','o')
% set(MapS,...
%     'DisplayName','S Cones',...
%     'MarkerFaceColor',[0 .45 .74],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',3,...
%     'Marker','o')
% 
% %Plot the weighted center and surround inputs for this particular cell
% figure;
% subplot(1,2,1);
% scatter(SortedTotalMosaic(1:ConesToBipolarCenter,1),SortedTotalMosaic(1:ConesToBipolarCenter,2),30,NormalizedCenterKernel,'fill')
% axis([-max(CenterX)*1.1 max(CenterX)*1.1 -max(CenterX)*1.1 max(CenterX)*1.1]);
% axis square
% title(strcat('Center (',num2str(ConesToBipolarCenter),' cones)'));
% subplot(1,2,2);
% scatter(SortedTotalMosaic(1:ConesToBipolarSurround,1),SortedTotalMosaic(1:ConesToBipolarSurround,2),30,NormalizedSurroundKernel,'fill')
% axis([-max(SurroundX)*1.1 max(SurroundX)*1.1 -max(SurroundX)*1.1 max(SurroundX)*1.1]);
% axis square
% title(strcat('Surround (',num2str(ConesToBipolarSurround),' cones)'));
% 
% %Plot the entire patch to see the receptive field in context (center,
% %surround, DoG)
% figure;
% subplot(2,3,1);
% scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,fullfieldCenter,'fill')
% axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
% axis square
% title('Center')
% subplot(2,3,2);
% scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,fullfieldSurround,'fill')
% axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
% axis square
% title('Surround')
% subplot(2,3,3)
% scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,DoG,'fill')
% axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
% axis square
% title('C-S')
% 
% %Plot the same thing, but visualize the weights of the cone classes
% subplot(2,3,4);
% hold on;
% CLPlot=plot3(LCoords(:,1),LCoords(:,2),CenterLCoordWeights,'r.');
% CMPlot=plot3(MCoords(:,1),MCoords(:,2),CenterMCoordWeights,'g.');
% CSPlot=plot3(SCoords(:,1),SCoords(:,2),CenterSCoordWeights,'b.');
% view([33 6]);
% title(strcat('Center (',num2str(cnL),'L/',num2str(cnM),'M/',num2str(cnS),'S)'));
% set(CLPlot,...
%     'DisplayName','Surround L',...
%     'MarkerFaceColor',[.64 .08 .18],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(CMPlot,...
%     'DisplayName','Surround M',...
%     'MarkerFaceColor',[.47 .67 .19],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(CSPlot,...
%     'DisplayName','Surround S',...
%     'MarkerFaceColor',[0 .45 .74],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% 
% subplot(2,3,5);
% hold on;
% SLPlot=plot3(LCoords(:,1),LCoords(:,2),SurroundLCoordWeights,'ro');
% SMPlot=plot3(MCoords(:,1),MCoords(:,2),SurroundMCoordWeights,'go');
% SSPlot=plot3(SCoords(:,1),SCoords(:,2),SurroundSCoordWeights,'bo');
% view([33 6]);
% title(strcat('Surround (',num2str(snL),'L/',num2str(snM),'M/',num2str(snS),'S)'));
% set(SLPlot,...
%     'DisplayName','Surround L',...
%     'MarkerFaceColor',[.64 .08 .18],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(SMPlot,...
%     'DisplayName','Surround M',...
%     'MarkerFaceColor',[.47 .67 .19],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(SSPlot,...
%     'DisplayName','Surround S',...
%     'MarkerFaceColor',[0 .45 .74],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% 
% subplot(2,3,6);
% hold on;
% DoGLPlot=plot3(LCoords(:,1),LCoords(:,2),LDoG,'r.');
% DoGMPlot=plot3(MCoords(:,1),MCoords(:,2),MDoG,'g.');
% DoGSPlot=plot3(SCoords(:,1),SCoords(:,2),SDoG,'b.');
% view([33 6]);
% set(DoGLPlot,...
%     'DisplayName','Surround L',...
%     'MarkerFaceColor',[.64 .08 .18],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(DoGMPlot,...
%     'DisplayName','Surround M',...
%     'MarkerFaceColor',[.47 .67 .19],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
% set(DoGSPlot,...
%     'DisplayName','Surround S',...
%     'MarkerFaceColor',[0 .45 .74],...
%     'LineStyle','none',...
%     'Color','none',...
%     'MarkerSize',2,...
%     'Marker','o')
