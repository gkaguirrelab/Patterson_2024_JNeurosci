%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateMidgets(Em)
%% GENERATE THE PATCH
%Load a universal hex grid
cd('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid');
load('HexGrid');

%Note the standard-issue distance between two points
ganglion_coords=pixel_coords;
PDist=pdist2(ganglion_coords(1,:),ganglion_coords(2,:));

%Jitter the cartesian coordinates to make a unique mosaic
ganglion_coords(:,1)=jitter(ganglion_coords(:,1));
ganglion_coords(:,2)=jitter(ganglion_coords(:,2));

%Determine the midget DF radius given the eccentricity
DFRadius=0.002738*Em.^1.327; %From 'cones per midget gc data' Excel file

%Set the distance between points (~2 x DFRadius(Em))
RFactor=(DFRadius*2)/PDist;
ganglion_coordsADJ=ganglion_coords*RFactor;
gDistances=pdist2(ganglion_coordsADJ,[0,0]);

%%
%Set up some parameters for the size and strength of the center and
%surround responses
CenterRadius=0.002738*Em.^1.327; %in mm
SurroundRadius=6*CenterRadius; %simplification of Kroner/Kaplan 1995
% SurroundRadius=0.7746*(CenterRadius.^0.3557); %Kroner/Kaplan Figure 4 fit
% SurroundRadius=11.67*(CenterRadius.^2)+(1.039)*CenterRadius+0.02969; %JOV 2002 fit
CSRadRatio=SurroundRadius/CenterRadius;
CSRadRatio=6;
CSStrengthRatio=.75;

% %What is the cone count for center and surround, given our particular location?
% ConesToCenter=floor(0.3736*Em^2 + 0.01126*Em + 1.026);
% ConesToSurround=(CSRadRatio^2)*ConesToCenter;
% 
% %Assign a sigma to determine the shape of the kernels
% CenterX=CenterRadius;
% SurroundX=SurroundRadius;
% cSigma=max(CenterX); %(we assume ~3 SDs make an RF radius)
% sSigma=max(SurroundX);
% 
% %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
% %from the origin
% CenterKernel=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
% NormalizedCenterKernel=CenterKernel/(sum(CenterKernel));
% SurroundKernel=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
% NormalizedSurroundKernel=(SurroundKernel/(sum(SurroundKernel)));
% 
% %Contextualize the kernels within the mosaic at large...
% fullfieldCenter=zeros(length(SortedTotalMosaic),1);
% fullfieldCenter(1:ConesToCenter,:)=fullfieldCenter(1:ConesToCenter,:)+NormalizedCenterKernel;
% fullfieldSurround=zeros(length(SortedTotalMosaic),1);
% fullfieldSurround(1:ConesToSurround,:)=fullfieldSurround(1:ConesToSurround,:)+NormalizedSurroundKernel;
% 
% %...and combine them to see the complete receptive field
% DoG=fullfieldCenter-(CSStrengthRatio*fullfieldSurround);
% 
% %Assign each kernel weight to an individual cone
% CenterLCoordWeights=MosaicL.*fullfieldCenter;
% cpL=sum(MosaicL.*fullfieldCenter);
% CenterMCoordWeights=MosaicM.*fullfieldCenter;
% cpM=sum(MosaicM.*fullfieldCenter);
% CenterSCoordWeights=MosaicS.*fullfieldCenter;
% cpS=sum(MosaicS.*fullfieldCenter);
% 
% SurroundLCoordWeights=MosaicL.*fullfieldSurround;
% spL=sum(MosaicL.*fullfieldSurround);
% SurroundMCoordWeights=MosaicM.*fullfieldSurround;
% spM=sum(MosaicM.*fullfieldSurround);
% SurroundSCoordWeights=MosaicS.*fullfieldSurround;
% spS=sum(MosaicS.*fullfieldSurround);
% 
% LDoG=CenterLCoordWeights-SurroundLCoordWeights;
% MDoG=CenterMCoordWeights-SurroundMCoordWeights;
% SDoG=CenterLCoordWeights-SurroundSCoordWeights;
% 
% %Keep track of the discrete cone counts
% cnL=sum(MosaicL(1:ConesToCenter));
% cnM=sum(MosaicM(1:ConesToCenter));
% cnS=sum(MosaicS(1:ConesToCenter));
% 
% snL=sum(MosaicL(1:ConesToSurround));
% snM=sum(MosaicM(1:ConesToSurround));
% snS=sum(MosaicS(1:ConesToSurround));
%%

%Capture cones in each ganglion cell DF
MidgetArraySize=50;
for g=1:MidgetArraySize
    Label=strcat('Cell',num2str(g));
    gX=ganglion_coordsADJ(g,1); %X coord of the ganglion cell
    gY=ganglion_coordsADJ(g,2); %Y coord of the ganglion cell
    MidgetCells.(matlab.lang.makeValidName(Label)).Coords=[gX,gY];
    
    %Populate the centers
    for f=1:length(SortedTotalMosaic(1:5000,:))
        if ((SortedTotalMosaic(f,1)-gX)/CenterRadius)^2+((SortedTotalMosaic(f,2)-gY)/CenterRadius)^2<=1
            CenterX(f,1)=SortedTotalMosaic(f,1);
            CenterX(f,2)=SortedTotalMosaic(f,2);
            CenterX(f,3)=pdist2([SortedTotalMosaic(f,1),SortedTotalMosaic(f,2)],[gX,gY]);
            CenterLCoords(f,1)=LCoords(f,1);
            CenterLCoords(f,2)=LCoords(f,2);
            CenterLCoords(f,3)=pdist2([LCoords(f,1),LCoords(f,2)],[gX,gY]);
            CenterMCoords(f,1)=MCoords(f,1);
            CenterMCoords(f,2)=MCoords(f,2);
            CenterMCoords(f,3)=pdist2([MCoords(f,1),MCoords(f,2)],[gX,gY]);
        end
    end

    CenterIndex=find(CenterX(:,1));
    CenterX=CenterX(CenterIndex,:);
    
    CenterLIndex=find(CenterLCoords(:,1));
    CenterLCoords=CenterLCoords(CenterLIndex,:);

    CenterMIndex=find(CenterMCoords(:,1));
    CenterMCoords=CenterMCoords(CenterMIndex,:);
    
    %Populate the surrounds
    for f=1:length(SortedTotalMosaic(1:5000,:))
        if ((SortedTotalMosaic(f,1)-gX)/SurroundRadius)^2+((SortedTotalMosaic(f,2)-gY)/SurroundRadius)^2<=1
            SurroundX(f,1)=SortedTotalMosaic(f,1);
            SurroundX(f,2)=SortedTotalMosaic(f,2);
            SurroundX(f,3)=pdist2([SortedTotalMosaic(f,1),SortedTotalMosaic(f,2)],[gX,gY]);
            SurroundLCoords(f,1)=LCoords(f,1);
            SurroundLCoords(f,2)=LCoords(f,2);
            SurroundLCoords(f,3)=pdist2([LCoords(f,1),LCoords(f,2)],[gX,gY]);
            SurroundMCoords(f,1)=MCoords(f,1);
            SurroundMCoords(f,2)=MCoords(f,2);
            SurroundMCoords(f,3)=pdist2([MCoords(f,1),MCoords(f,2)],[gX,gY]);
        end
    end
    
    SurroundIndex=find(SurroundX(:,1));
    SurroundX=SurroundX(SurroundIndex,:);
    
    SurroundLIndex=find(SurroundLCoords(:,1));
    SurroundLCoords=SurroundLCoords(SurroundLIndex,:);

    SurroundMIndex=find(SurroundMCoords(:,1));
    SurroundMCoords=SurroundMCoords(SurroundMIndex,:);
    
    %Assign a center sigma to determine the shape of the kernels
    cSigma=max(CenterX(:,3)); %(we assume ~1 SD makes an RF radius)
    sSigma=max(SurroundX(:,3));

    %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
    %from the origin
    CenterKernel=exp((-1/(2*cSigma^2))*(abs(CenterX(:,3)).^2));
    NormalizedCenterKernel=CenterKernel/(sum(CenterKernel));
    SurroundKernel=exp((-1/(2*sSigma^2))*(abs(SurroundX(:,3)).^2));
    NormalizedSurroundKernel=(SurroundKernel/(sum(SurroundKernel)));
    
    MidgetCells.(matlab.lang.makeValidName(Label)).Center.CenterConeCoords=CenterX;
    MidgetCells.(matlab.lang.makeValidName(Label)).Center.ConesToCenter=length(MidgetCells.(matlab.lang.makeValidName(Label)).Center.CenterConeCoords);
    MidgetCells.(matlab.lang.makeValidName(Label)).Center.NormKernel=NormalizedCenterKernel;
    
    MidgetCells.(matlab.lang.makeValidName(Label)).Center.LCones.Coords=CenterLCoords;
    MidgetCells.(matlab.lang.makeValidName(Label)).Center.MCones.Coords=CenterMCoords;
    
    
    MidgetCells.(matlab.lang.makeValidName(Label)).Surround.ConeCoords.All=SurroundX;
    MidgetCells.(matlab.lang.makeValidName(Label)).Surround.ConeCoords.LCones=SurroundLCoords;
    MidgetCells.(matlab.lang.makeValidName(Label)).Surround.ConeCoords.MCones=SurroundMCoords;
    MidgetCells.(matlab.lang.makeValidName(Label)).Surround.ConesToSurround=length(MidgetCells.(matlab.lang.makeValidName(Label)).Surround.ConeCoords.All);
    MidgetCells.(matlab.lang.makeValidName(Label)).Surround.NormKernel=NormalizedSurroundKernel;
end

%%

figure;
hold on;
MapL=plot(LCoords(:,1),LCoords(:,2),'r.');
MapM=plot(MCoords(:,1),MCoords(:,2),'g.');
MapS=plot(SCoords(:,1),SCoords(:,2),'b.');
axis square
set(gca, 'color', [0.25 0.25 0.25]);
set(MapL,...
    'DisplayName','L Cones',...
    'MarkerFaceColor',[.64 .08 .18],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')
set(MapM,...
    'DisplayName','M Cones',...
    'MarkerFaceColor',[.47 .67 .19],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')
set(MapS,...
    'DisplayName','S Cones',...
    'MarkerFaceColor',[0 .45 .74],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')

for g=1:MidgetArraySize
    [X,Y]=circle([ganglion_coordsADJ(g,1),ganglion_coordsADJ(g,2)],DFRadius,1000);
    plot(X,Y,'b');
    text(ganglion_coordsADJ(g,1),ganglion_coordsADJ(g,2),num2str(g),'Color','white');
end