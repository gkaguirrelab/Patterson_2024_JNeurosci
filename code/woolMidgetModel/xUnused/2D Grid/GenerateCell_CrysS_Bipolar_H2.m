%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateCell_CrysS_Bipolar_H2(Em)
%Assign the cell as an ON or an OFF midget cell
CellType=((randi(2)-2)*2)+1;

%% GENERATE THE CONE MOSAIC & BIPOLAR LAYER
[ConeMosaic,BipolarLayer]=GenerateBipolarLayer_CrysS_ONvOFF(Em,CellType);

%% GENERATE THE CELL

%What is the cone count for center and surround, given our particular location?
CSStrengthRatio=1; %How much weaker is the surround than the center?
ConesToCenter=round(0.3736*Em^2 + 0.01126*Em + 1.026);

%Keep track of the discrete cone counts
cnL=sum(~isnan(ConeMosaic.LConeCoords(1:ConesToCenter)));
cnM=sum(~isnan(ConeMosaic.MConeCoords(1:ConesToCenter)));
cnS=sum(~isnan(ConeMosaic.SConeCoords(1:ConesToCenter)));

SurroundIndex=unique(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter));
snL=sum(logical(~isnan(ConeMosaic.LConeCoords(unique(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter))))));
snM=sum(logical(~isnan(ConeMosaic.MConeCoords(unique(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter))))));
snS=sum(logical(~isnan(ConeMosaic.SConeCoords(unique(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter))))));

ConesToSurround=snL+snM+snS;

%Assign the center and surround cone coordinates, starting from (0,0) and going 
%outward until the cone count is fulfilled (this ensures that the center
%and surround share the centralmost cones)
CenterX=ConeMosaic.AllConeDistances(1:ConesToCenter);
SurroundX=ConeMosaic.AllConeDistances(SurroundIndex);

%Assign a sigma to determine the shape of the kernels
cSigma=max(CenterX); %(we assume ~1 SD make an RF radius)
sSigma=max(SurroundX);

%Apply a normalized rbfGaussian kernel to each cone as a function of its distance
%from the origin
CenterKernel=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
NormalizedCenterKernel=CenterKernel/(sum(CenterKernel));

%Contextualize the kernels within the mosaic at large...
fullfieldCenter=zeros(length(ConeMosaic.AllConeCoords),1);
fullfieldCenter(1:ConesToCenter,:)=fullfieldCenter(1:ConesToCenter,:)+NormalizedCenterKernel;
fullfieldSurround=zeros(length(ConeMosaic.AllConeCoords),1);

LIndex=logical(~isnan(ConeMosaic.LConeCoords(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter)))).*BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter);
MIndex=logical(~isnan(ConeMosaic.MConeCoords(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter)))).*BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter);
SIndex=logical(~isnan(ConeMosaic.SConeCoords(BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter)))).*BipolarLayer.SurroundArray.Indices(:,1:ConesToCenter);

%If it's an ON cell, zero out any S cone contribution to the center
%(the surround has already been taken care of during bipolar generation)
if CellType==1
    fullfieldCenter(find(SIndex(1,:)))=0;
end

for g=1:ConesToCenter
    AddL=LIndex(:,g);
    AddL(AddL==0)=[];
    LWeights=BipolarLayer.SurroundArray.LConeWeights(:,g)*fullfieldCenter(g);
    LWeights(LWeights==0)=[];
    if isempty(LWeights)
        LWeights=0;
    end
    
    AddM=MIndex(:,g);
    AddM(AddM==0)=[];
    MWeights=BipolarLayer.SurroundArray.MConeWeights(:,g)*fullfieldCenter(g);
    MWeights(MWeights==0)=[];
    if isempty(MWeights)
        MWeights=0;
    end
    
    AddS=SIndex(:,g);
    AddS(AddS==0)=[];
    SWeights=BipolarLayer.SurroundArray.SConeWeights(:,g)*fullfieldCenter(g);
    SWeights(SWeights==0)=[];
    if isempty(SWeights)
        SWeights=0;
    end
    
    fullfieldSurround(AddL)=fullfieldSurround(AddL)+LWeights;
    fullfieldSurround(AddM)=fullfieldSurround(AddM)+MWeights;
    fullfieldSurround(AddS)=fullfieldSurround(AddS)+SWeights;
end

%...and combine them to see the complete receptive field
DoG=fullfieldCenter-(fullfieldSurround);

%Assign each kernel weight to an individual cone
CenterLCoordWeights=(~isnan(ConeMosaic.LConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter);
cpL=sum((~isnan(ConeMosaic.LConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter));
CenterMCoordWeights=(~isnan(ConeMosaic.MConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter);
cpM=sum((~isnan(ConeMosaic.MConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter));
CenterSCoordWeights=(~isnan(ConeMosaic.SConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter);
cpS=sum((~isnan(ConeMosaic.SConeCoords(1:ConesToCenter)))'.*fullfieldCenter(1:ConesToCenter));

%Total surround:
SurroundLCoordWeights=(~isnan(ConeMosaic.LConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex);
spL=sum((~isnan(ConeMosaic.LConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex));
SurroundMCoordWeights=(~isnan(ConeMosaic.MConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex);
spM=sum((~isnan(ConeMosaic.MConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex));
SurroundSCoordWeights=(~isnan(ConeMosaic.SConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex);
spS=sum((~isnan(ConeMosaic.SConeCoords(SurroundIndex))).*fullfieldSurround(SurroundIndex));


%Compute the DoGs
LDoG=CellType*padarray(CenterLCoordWeights,ConeMosaic.ConesPerMM-ConesToCenter,0,'post')-padarray(SurroundLCoordWeights,ConeMosaic.ConesPerMM-ConesToSurround,0,'post');
MDoG=CellType*padarray(CenterMCoordWeights,ConeMosaic.ConesPerMM-ConesToCenter,0,'post')-padarray(SurroundMCoordWeights,ConeMosaic.ConesPerMM-ConesToSurround,0,'post');
SDoG=CellType*padarray(CenterSCoordWeights,ConeMosaic.ConesPerMM-ConesToCenter,0,'post')-padarray(SurroundSCoordWeights,ConeMosaic.ConesPerMM-ConesToSurround,0,'post');

%% SAVE TO OUTPUT

% GLOBAL
Cell.RetinaWeights=ConeMosaic.RetinaWeights;
Cell.CenterConeCount=[cnL cnM cnS];
Cell.SurroundConeCount=[snL snM snS];
Cell.CenterWeights=[cpL cpM cpS];
Cell.SurroundWeights=[spL spM spS];
Cell.CSStrengthRatio=CSStrengthRatio;
Cell.Polarity=CellType;

Cell.Global.ConesPerMM=ConeMosaic.ConesPerMM;
Cell.Global.AllConeCoords=ConeMosaic.AllConeCoords;
Cell.Global.AllConeDistances=ConeMosaic.AllConeDistances;
Cell.Global.InterconeDistance=ConeMosaic.InterconeDistance;

Cell.Global.LConeCoords=ConeMosaic.LConeCoords;
Cell.Global.MConeCoords=ConeMosaic.MConeCoords;
Cell.Global.SConeCoords=ConeMosaic.SConeCoords;

Cell.Global.DoG.ReceptiveFieldKernel=DoG;
Cell.Global.DoG.LDoG=LDoG;
Cell.Global.DoG.MDoG=MDoG;
Cell.Global.DoG.SDoG=SDoG;

% CENTER
Cell.Center.ConesToCenter=ConesToCenter;
Cell.Center.CenterConeCoords=ConeMosaic.AllConeCoords(1:ConesToCenter,:);

Cell.Center.LCones.Count=cnL;
Cell.Center.LCones.Coords=ConeMosaic.LConeCoords(1:ConesToCenter,:);
Cell.Center.LCones.EachWeight=CenterLCoordWeights;
Cell.Center.LCones.TotalWeight=cpL;

Cell.Center.MCones.Count=cnM;
Cell.Center.MCones.Coords=ConeMosaic.MConeCoords(1:ConesToCenter,:);
Cell.Center.MCones.EachWeight=CenterMCoordWeights;
Cell.Center.MCones.TotalWeight=cpM;

Cell.Center.SCones.Count=cnS;
Cell.Center.SCones.Coords=ConeMosaic.SConeCoords(1:ConesToCenter,:);
Cell.Center.SCones.EachWeight=CenterSCoordWeights;
Cell.Center.SCones.TotalWeight=cpS;

Cell.Center.CenterConeDistances=CenterX;
Cell.Center.CenterSigma=cSigma;
Cell.Center.CenterKernel=fullfieldCenter;

% SURROUND
Cell.Surround.ConesToSurround=ConesToSurround;
Cell.Surround.SurroundConeCoords=ConeMosaic.AllConeCoords(SurroundIndex,:);

Cell.Surround.LCones.Count=snL;
Cell.Surround.LCones.Coords=ConeMosaic.LConeCoords(SurroundIndex,:);
Cell.Surround.LCones.EachWeight=SurroundLCoordWeights;
Cell.Surround.LCones.TotalWeight=spL;

Cell.Surround.MCones.Count=snM;
Cell.Surround.MCones.Coords=ConeMosaic.MConeCoords(SurroundIndex,:);
Cell.Surround.MCones.EachWeight=SurroundMCoordWeights;
Cell.Surround.MCones.TotalWeight=spM;

Cell.Surround.SCones.Count=snS;
Cell.Surround.SCones.Coords=ConeMosaic.SConeCoords(SurroundIndex,:);
Cell.Surround.SCones.EachWeight=SurroundSCoordWeights;
Cell.Surround.SCones.TotalWeight=spS;

Cell.Surround.SurroundConeDistances=SurroundX;
Cell.Surround.SurroundSigma=sSigma;
Cell.Surround.SurroundKernel=fullfieldSurround;

%% PLOT THE DATA
% 
% Plot the 1-mm2 patch and its L,M,S cones
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
% Plot the weighted center and surround inputs for this particular cell
% figure;
% subplot(1,2,1);
% scatter(SortedTotalMosaic(1:ConesToCenter,1),SortedTotalMosaic(1:ConesToCenter,2),30,NormalizedCenterKernel,'fill')
% axis([-max(CenterX)*1.1 max(CenterX)*1.1 -max(CenterX)*1.1 max(CenterX)*1.1]);
% axis square
% title(strcat('Center (',num2str(ConesToCenter),' cones)'));
% subplot(1,2,2);
% scatter(SortedTotalMosaic(1:ConesToSurround,1),SortedTotalMosaic(1:ConesToSurround,2),30,NormalizedSurroundKernel,'fill')
% axis([-max(SurroundX)*1.1 max(SurroundX)*1.1 -max(SurroundX)*1.1 max(SurroundX)*1.1]);
% axis square
% title(strcat('Surround (',num2str(ConesToSurround),' cones)'));
% 
% Plot the entire patch to see the receptive field in context (center,
% surround, DoG)
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
% Plot the same thing, but visualize the weights of the cone classes
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
