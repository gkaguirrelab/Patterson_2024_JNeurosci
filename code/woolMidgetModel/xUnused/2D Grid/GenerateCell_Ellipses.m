%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateCell_Ellipses(Em)
%% GENERATE THE PATCH
%Load a universal hex grid
cd('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid');
load('HexGrid');

%Jitter the cartesian coordinates to make a unique mosaic
pixel_coords(:,1)=jitter(pixel_coords(:,1));
pixel_coords(:,2)=jitter(pixel_coords(:,2));

%Determine the cone density given the eccentricity
ConesPerMM=ceil((52170*exp(-1.057.*Em))+(15550*exp(-0.1605.*Em))); %From 'cones per midget gc data' Excel file
DFRadius=0.002738*Em.^1.327; %From 'cones per midget gc data' Excel file

%Determine the aperture width needed to encompass the correct number of
%cones
r=1;
ConeCount=0;
while ConeCount<ConesPerMM
    r=r+.1;
    counter=0;
    inc_pix_coords=[];
    for p=1:length(pixel_coords)
        if pixel_coords(p,1)>-r && pixel_coords(p,1)<r
            if pixel_coords(p,2)>-r && pixel_coords(p,2)<r
                inc_pix_coords(end+1,1)=pixel_coords(p,1);
                inc_pix_coords(end,2)=pixel_coords(p,2);
                counter=counter+1;
            end
        end
    end
    ConeCount=counter;
end

%Adjust the aperture window to a 1-mm2 portion...now you have cones/mm2
RFactor=1/(2*r);
pixel_coordsADJ=pixel_coords*RFactor;
inc_pix_coordsADJ=inc_pix_coords*RFactor;

%Compute the distances of each cone from the center (0,0) and sort the
%coordinates from closest to farthest
distances=pdist2(inc_pix_coordsADJ,[0,0]);
SortedTotalMosaic=sortrows([inc_pix_coordsADJ,distances],3);

%Determine each cone's nearest neighbors and use their average as an estimate of
%intercone distance across the retina
NearestNeighbors=zeros(length(SortedTotalMosaic),1);
HowManyNeighbors=1;
for q=1:length(SortedTotalMosaic)
    dist=sort(pdist2(SortedTotalMosaic(q,:),SortedTotalMosaic));
    NearestNeighbors(q)=mean(dist(2:HowManyNeighbors+1));
end
InterConeDistance=mean(NearestNeighbors);

% Assign a retinal L:M ratio for this patch
% (L:M:S ratios vary across retinas, lognormally)
LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
RandLM=exp(random(LogLMDistribution));
pL=RandLM/(RandLM+1); %Likelihood of L
pM=1-pL; %Likelihood of M
pS=0.0; %Likelihood of S

if pS~=0
    pL=pL*(1-pS);
    pM=pM*(1-pS);
end

%Assign the patch to L, M, and S 'zones' (for assigning discrete cones later)
LLIndex=0;
LUIndex=floor(ConeCount*pL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+floor(ConeCount*pM);
SLIndex=MUIndex+1;
SUIndex=ConeCount;

%Generate the random assignments to the entire patch
MosaicAssign(:,1)=randi(ConeCount,1,length(SortedTotalMosaic));

%Determine which center values are Ls, Ms, and Ss...
MosaicL=MosaicAssign<=LUIndex;
MosaicM=(MosaicAssign>=MLIndex & MosaicAssign<=MUIndex);
MosaicS=(MosaicAssign>=SLIndex & MosaicAssign<=SUIndex);

%Give the cone assignments their corresponding coordinates
LCoords=[SortedTotalMosaic(:,1).*MosaicL SortedTotalMosaic(:,2).*MosaicL];
LCoords(LCoords==0)=NaN;
MCoords=[SortedTotalMosaic(:,1).*MosaicM SortedTotalMosaic(:,2).*MosaicM];
MCoords(MCoords==0)=NaN;
SCoords=[SortedTotalMosaic(:,1).*MosaicS SortedTotalMosaic(:,2).*MosaicS];
SCoords(SCoords==0)=NaN;

% %Quasicrystalline S addition: Use this if you are not going to depend on
% %random generation of S locations as written above (be sure pS=0)
% pS=.1;
% SDistance=InterConeDistance*2.5;
% 
% NearestNeighbors=zeros(length(SortedTotalMosaic));
% for q=1:length(SortedTotalMosaic)
% dist=pdist2(SortedTotalMosaic(q,:),SortedTotalMosaic);
% NearestNeighbors(q,:)=dist;
% end
% 
% done=0;
% qualify=[];
% qualify(1)=randi(10);
% vector=[];
% vector=find(NearestNeighbors(qualify(1):end,qualify(1))>SDistance);
% qualify(2)=vector(1)+qualify(1)-1;
% 
% for x=2:length(SortedTotalMosaic)
%     if x>length(qualify)
%         break
%     end
%     vector=find(NearestNeighbors(qualify(x):end,qualify(x))>SDistance);
%     for y=1:length(vector)
%         if NearestNeighbors(vector(y)+qualify(x)-1,qualify(1:x-1))>SDistance
%             qualify(x+1)=vector(y)+qualify(x)-1;
%             vector=[];
%             break
%         end
%     end
% end
% 
% CrystallineSCones=SortedTotalMosaic(qualify',1:2);
% for g=1:length(SortedTotalMosaic)
%     z(g,:)=ismember(SortedTotalMosaic(g,1:2),CrystallineSCones);
% end
% MosaicS=z(:,1);
% SCoords=[SortedTotalMosaic(:,1).*MosaicS SortedTotalMosaic(:,2).*MosaicS];
% SCoords(SCoords==0)=NaN;

%% GENERATE THE CELL

%Set up some parameters for the size and strength of the center and
%surround responses
CenterRadius=0.002738*Em.^1.327; %in mm
SurroundRadius=6*CenterRadius; %simplification of Kroner/Kaplan 1995
% SurroundRadius=0.7746*(CenterRadius.^0.3557); %Kroner/Kaplan Figure 4 fit
% SurroundRadius=11.67*(CenterRadius.^2)+(1.039)*CenterRadius+0.02969; %JOV 2002 fit
CSRadRatio=SurroundRadius/CenterRadius;
CSStrengthRatio=.75;

%What is the cone count for center and surround, given our particular location?
ConesToCenter=ceil(0.3736*Em^2 + 0.01126*Em + 1.026);
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Build a ~2:1 ellipse to capture cones
CenterEllipseRadiusA=0.003835*Em.^1.327; %in mm
CenterEllipseRadiusB=0.001952*Em.^1.328; %in mm
[CEllipseX,CEllipseY]=ellipse(CenterEllipseRadiusA,CenterEllipseRadiusB,0,0,0);

for f=1:length(SortedTotalMosaic(1:50,:))
    if (SortedTotalMosaic(f,1)/CenterEllipseRadiusA)^2+(SortedTotalMosaic(f,2)/CenterEllipseRadiusB)^2<=1
        CenterX(f,1)=SortedTotalMosaic(f,1);
        CenterX(f,2)=SortedTotalMosaic(f,2);
        CenterX(f,3)=distances(f,:);
    end
end
CEllipseIndex=find(CenterX(:,1));
CenterX=CenterX(CEllipseIndex,:);
ConesToCenter=length(CenterX);

SurroundEllipseRadiusA=CenterEllipseRadiusA*CSRadRatio; %in mm
SurroundEllipseRadiusB=CenterEllipseRadiusB*CSRadRatio; %in mm
[SEllipseX,SEllipseY]=ellipse(SurroundEllipseRadiusA,SurroundEllipseRadiusB,0,0,0);

for f=1:length(SortedTotalMosaic(:,1))
    if (SortedTotalMosaic(f,1)/SurroundEllipseRadiusA)^2+(SortedTotalMosaic(f,2)/SurroundEllipseRadiusB)^2<=1
        SurroundX(f,1)=SortedTotalMosaic(f,1);
        SurroundX(f,2)=SortedTotalMosaic(f,2);
        SurroundX(f,3)=distances(f,:);
    end
end
SEllipseIndex=find(SurroundX(:,1));
SurroundX=SurroundX(SEllipseIndex,:);
ConesToSurround=length(SurroundX);

%Assign a sigma to determine the shape of the kernels
cSigmaA=CenterEllipseRadiusA; %(we assume ~1 SD make an RF radius)
cSigmaB=CenterEllipseRadiusB;
sSigmaA=SurroundEllipseRadiusA; %(we assume ~1 SD make an RF radius)
sSigmaB=SurroundEllipseRadiusB;

%Apply a normalized rbfGaussian kernel to each cone as a function of its distance
%from the origin
CenterEllipticalKernel=exp(-(((CenterX(:,1))./cSigmaA).^2+((CenterX(:,2))./cSigmaB).^2));
NormalizedCenterEllipticalKernel=CenterEllipticalKernel/(sum(CenterEllipticalKernel));
SurroundEllipticalKernel=exp(-(((SurroundX(:,1))./sSigmaA).^2+((SurroundX(:,2))./sSigmaB).^2));
NormalizedSurroundEllipticalKernel=SurroundEllipticalKernel/(sum(SurroundEllipticalKernel));

%Contextualize the kernels within the mosaic at large...
fullfieldCenter=zeros(length(SortedTotalMosaic),1);
fullfieldCenter(CEllipseIndex,:)=fullfieldCenter(CEllipseIndex,:)+NormalizedCenterEllipticalKernel;
fullfieldSurround=zeros(length(SortedTotalMosaic),1);
fullfieldSurround(SEllipseIndex,:)=fullfieldSurround(SEllipseIndex,:)+NormalizedSurroundEllipticalKernel;

%...and combine them to see the complete receptive field
DoG=fullfieldCenter-(CSStrengthRatio*fullfieldSurround);

%Assign each kernel weight to an individual cone
CenterLCoordWeights=MosaicL.*fullfieldCenter;
cpL=sum(MosaicL.*fullfieldCenter);
CenterMCoordWeights=MosaicM.*fullfieldCenter;
cpM=sum(MosaicM.*fullfieldCenter);
CenterSCoordWeights=MosaicS.*fullfieldCenter;
cpS=sum(MosaicS.*fullfieldCenter);

SurroundLCoordWeights=MosaicL.*fullfieldSurround;
spL=sum(MosaicL.*fullfieldSurround);
SurroundMCoordWeights=MosaicM.*fullfieldSurround;
spM=sum(MosaicM.*fullfieldSurround);
SurroundSCoordWeights=MosaicS.*fullfieldSurround;
spS=sum(MosaicS.*fullfieldSurround);

LDoG=CenterLCoordWeights-SurroundLCoordWeights;
MDoG=CenterMCoordWeights-SurroundMCoordWeights;
SDoG=CenterLCoordWeights-SurroundSCoordWeights;

%Keep track of the discrete cone counts
cnL=sum(MosaicL(1:ConesToCenter));
cnM=sum(MosaicM(1:ConesToCenter));
cnS=sum(MosaicS(1:ConesToCenter));

snL=sum(MosaicL(1:ConesToSurround));
snM=sum(MosaicM(1:ConesToSurround));
snS=sum(MosaicS(1:ConesToSurround));

%% SAVE TO OUTPUT

% GLOBAL
Cell.RetinaWeights=[pL pM pS];
Cell.CenterConeCount=[cnL cnM cnS];
Cell.SurroundConeCount=[snL snM snS];
Cell.CenterWeights=[cpL cpM cpS];
Cell.SurroundWeights=[spL spM spS];

Cell.Global.ConesPerMM=ConeCount;
Cell.Global.AllConeCoords=SortedTotalMosaic;
Cell.Global.AllConeDistances=distances;
Cell.Global.InterconeDistance=InterConeDistance;

Cell.Global.LConeCoords=LCoords;
Cell.Global.MConeCoords=MCoords;
Cell.Global.SConeCoords=SCoords;

Cell.Global.DoG.ReceptiveFieldKernel=DoG;
Cell.Global.DoG.LDoG=LDoG;
Cell.Global.DoG.MDoG=MDoG;
Cell.Global.DoG.SDoG=SDoG;

% CENTER
Cell.Center.ConesToCenter=ConesToCenter;
Cell.Center.CenterConeCoords=SortedTotalMosaic(CEllipseIndex,:);
Cell.Center.EllipseIndex=CEllipseIndex;

Cell.Center.LCones.Count=cnL;
Cell.Center.LCones.Coords=LCoords(CEllipseIndex,:);
Cell.Center.LCones.EachWeight=CenterLCoordWeights;
Cell.Center.LCones.TotalWeight=cpL;

Cell.Center.MCones.Count=cnM;
Cell.Center.MCones.Coords=MCoords(CEllipseIndex,:);
Cell.Center.MCones.EachWeight=CenterMCoordWeights;
Cell.Center.MCones.TotalWeight=cpM;

Cell.Center.SCones.Count=cnS;
Cell.Center.SCones.Coords=SCoords(CEllipseIndex,:);
Cell.Center.SCones.EachWeight=CenterSCoordWeights;
Cell.Center.SCones.TotalWeight=cpS;

Cell.Center.CenterConeDistances=CenterX(:,3);
Cell.Center.CenterSigmas=[cSigmaA cSigmaB];
Cell.Center.CenterKernel=NormalizedCenterEllipticalKernel;
Cell.Center.CenterKernelInField=fullfieldCenter;

% SURROUND
Cell.Surround.ConesToSurround=ConesToSurround;
Cell.Surround.SurroundConeCoords=SortedTotalMosaic(SEllipseIndex,:);
Cell.Surround.EllipseIndex=SEllipseIndex;

Cell.Surround.LCones.Count=snL;
Cell.Surround.LCones.Coords=LCoords(SEllipseIndex,:);
Cell.Surround.LCones.EachWeight=SurroundLCoordWeights;
Cell.Surround.LCones.TotalWeight=spL;

Cell.Surround.MCones.Count=snM;
Cell.Surround.MCones.Coords=MCoords(SEllipseIndex,:);
Cell.Surround.MCones.EachWeight=SurroundMCoordWeights;
Cell.Surround.MCones.TotalWeight=spM;

Cell.Surround.SCones.Count=snS;
Cell.Surround.SCones.Coords=SCoords(SEllipseIndex,:);
Cell.Surround.SCones.EachWeight=SurroundSCoordWeights;
Cell.Surround.SCones.TotalWeight=spS;

Cell.Surround.SurroundConeDistances=SurroundX(:,3);
Cell.Surround.SurroundSigmas=[sSigmaA sSigmaB];
Cell.Surround.SurroundKernel=NormalizedSurroundEllipticalKernel;
Cell.Surround.SurroundKernelInField=fullfieldSurround;

%% PLOT THE DATA

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
