%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateCell_CrysS_H2(Em)
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
parfor q=1:length(SortedTotalMosaic)
    dist=sort(pdist2(SortedTotalMosaic(q,:),SortedTotalMosaic));
    NearestNeighbors(q)=mean(dist(2:HowManyNeighbors+1));
end
InterConeDistance=mean(NearestNeighbors);

% Assign a retinal L:M ratio for this patch
% (L:M:S ratios vary across retinas, lognormally)
LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
RandLM=exp(random(LogLMDistribution));
rpL=RandLM/(RandLM+1); %Likelihood of L
rpM=1-rpL; %Likelihood of M
rpS=0; %Likelihood of S

if rpS~=0
    rpL=rpL*(1-rpS);
    rpM=rpM*(1-rpS);
end

%Assign the patch to L, M, and S 'zones' (for assigning discrete cones later)
LLIndex=0;
LUIndex=floor(ConeCount*rpL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+floor(ConeCount*rpM);
SLIndex=MUIndex+1;
SUIndex=ConeCount;

%Generate the random assignments to the entire patch
MosaicAssign(:,1)=randi(ConeCount,1,length(SortedTotalMosaic));

%Determine which center values are Ls, Ms, and Ss...
MosaicL=MosaicAssign<=LUIndex;
MosaicM=(MosaicAssign>=MLIndex & MosaicAssign<=MUIndex);
MosaicS=(MosaicAssign>=SLIndex & MosaicAssign<=SUIndex);

%Quasicrystalline S addition: Use this if you are not going to depend on
%random generation of S locations as written above (be sure pS=0)
rpS=0.1;
SDistance=InterConeDistance*sqrt(100/(rpS*100));

%Compute the distance of all cones to each cone
NearestNeighbors=zeros(length(SortedTotalMosaic));
parfor q=1:length(SortedTotalMosaic)
dist=pdist2(SortedTotalMosaic(q,:),SortedTotalMosaic);
NearestNeighbors(q,:)=dist;
end

%Find the cones that are all at least SDistance apart
done=0;
qualify=[];
qualify(1)=randi(10);
vector=[];
vector=find(NearestNeighbors(qualify(1):end,qualify(1))>SDistance);
qualify(2)=vector(1)+qualify(1)-1;

for x=2:length(SortedTotalMosaic)
    if x>length(qualify)
        break
    end
    vector=find(NearestNeighbors(qualify(x):end,qualify(x))>SDistance);
    for y=1:length(vector)
        if NearestNeighbors(vector(y)+qualify(x)-1,qualify(1:x-1))>SDistance
            qualify(x+1)=vector(y)+qualify(x)-1;
            vector=[];
            break
        end
    end
end

%Log the cone coordinates and create a logical S locations vector
CrystallineSCones=SortedTotalMosaic(qualify',1:2);
parfor g=1:length(SortedTotalMosaic)
    z(g,:)=ismember(SortedTotalMosaic(g,1:2),CrystallineSCones);
end
MosaicS=z(:,1);
SCoords=[SortedTotalMosaic(:,1).*MosaicS SortedTotalMosaic(:,2).*MosaicS];
SCoords(SCoords==0)=NaN;

%Adjust the L and M vectors since you took some of their elements over for
%S cones
SOverlapL=MosaicL-MosaicS;
SOverlapL(SOverlapL==-1)=0;
MosaicL=SOverlapL;

SOverlapM=MosaicM-MosaicS;
SOverlapM(SOverlapM==-1)=0;
MosaicM=SOverlapM;

%Give the cone assignments their corresponding coordinates
LCoords=[SortedTotalMosaic(:,1).*MosaicL SortedTotalMosaic(:,2).*MosaicL];
LCoords(LCoords==0)=NaN;
MCoords=[SortedTotalMosaic(:,1).*MosaicM SortedTotalMosaic(:,2).*MosaicM];
MCoords(MCoords==0)=NaN;
SCoords=[SortedTotalMosaic(:,1).*MosaicS SortedTotalMosaic(:,2).*MosaicS];
SCoords(SCoords==0)=NaN;

%Report the new effective pL, pM, pS for this retina
pL=sum(MosaicL)/length(SortedTotalMosaic);
pM=sum(MosaicM)/length(SortedTotalMosaic);
pS=sum(MosaicS)/length(SortedTotalMosaic);

%% GENERATE THE CELL

%What is the cone count for center and surround, given our particular location?
CSStrengthRatio=1; %How much weaker is the surround than the center?
CSRadRatio=6; %Assume that H1 and H2 cells have similar breadth (i.e., similar surround dimensions)
ConesToCenter=round(0.3736*Em^2 + 0.01126*Em + 1.026);
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Keep track of the discrete cone counts
cnL=sum(MosaicL(1:ConesToCenter));
cnM=sum(MosaicM(1:ConesToCenter));
cnS=sum(MosaicS(1:ConesToCenter));

snL=sum(MosaicL(1:ConesToSurround));
snM=sum(MosaicM(1:ConesToSurround));
snS=sum(MosaicS(1:ConesToSurround));

%Set up some parameters for the size and strength of the center and
%surround responses

H1StrengthRatio=.825; %H1 Cells make up ~82.5% of all horizontal synapses
if cnL+cnM==0
    H1LMWeight=0; %If there are no L/M cones in the RF center, then there is no H1 contribution to the RF
else
    H1LMWeight=1; %H1 cells only connect to L and M cones...
end
H1SWeight=0; % ...and never S cones :(

H2StrengthRatio=.175; %H2 cells make up 17.5% of all horizontal synapses
H2LMWeight=1/6; %H2 signaling from a single L/M cone is 1/6 the strength of that from a single S cone
H2SWeight=1;


%Assign the center and surround cone coordinates, starting from (0,0) and going 
%outward until the cone count is fulfilled (this ensures that the center
%and surround share the centralmost cones)
CenterX=distances(1:ConesToCenter);
SurroundX=distances(1:ConesToSurround);

%Assign a sigma to determine the shape of the kernels
cSigma=max(CenterX); %(we assume ~1 SD make an RF radius)
sSigma=max(SurroundX);

%Apply a normalized rbfGaussian kernel to each cone as a function of its distance
%from the origin
CenterKernel=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
NormalizedCenterKernel=CenterKernel/(sum(CenterKernel));
SurroundKernel=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
NormalizedSurroundKernel=(SurroundKernel/(sum(SurroundKernel)));

%Contextualize the kernels within the mosaic at large...
fullfieldCenter=zeros(length(SortedTotalMosaic),1);
fullfieldCenter(1:ConesToCenter,:)=fullfieldCenter(1:ConesToCenter,:)+NormalizedCenterKernel;
fullfieldSurround=zeros(length(SortedTotalMosaic),1);
fullfieldSurround(1:ConesToSurround,:)=fullfieldSurround(1:ConesToSurround,:)+NormalizedSurroundKernel;

%...and combine them to see the complete receptive field
DoG=fullfieldCenter-(H1StrengthRatio*fullfieldSurround);

%Assign each kernel weight to an individual cone
CenterLCoordWeights=MosaicL.*fullfieldCenter;
cpL=sum(MosaicL.*fullfieldCenter);
CenterMCoordWeights=MosaicM.*fullfieldCenter;
cpM=sum(MosaicM.*fullfieldCenter);
CenterSCoordWeights=MosaicS.*fullfieldCenter;
cpS=sum(MosaicS.*fullfieldCenter);

%For the H1 surround:
H1SurroundLCoordWeights=MosaicL.*fullfieldSurround*H1LMWeight*H1StrengthRatio;
H1SurroundMCoordWeights=MosaicM.*fullfieldSurround*H1LMWeight*H1StrengthRatio;
H1SurroundSCoordWeights=MosaicS.*fullfieldSurround*H1SWeight*H1StrengthRatio;

%For the H2 surround:
H2SurroundLCoordWeights=MosaicL.*fullfieldSurround*H2LMWeight*H2StrengthRatio;
H2SurroundMCoordWeights=MosaicM.*fullfieldSurround*H2LMWeight*H2StrengthRatio;
H2SurroundSCoordWeights=MosaicS.*fullfieldSurround*H2SWeight*H2StrengthRatio;

%Total surround:
SurroundLCoordWeights=H1SurroundLCoordWeights+H2SurroundLCoordWeights;
SurroundMCoordWeights=H1SurroundMCoordWeights+H2SurroundMCoordWeights;
SurroundSCoordWeights=H1SurroundSCoordWeights+H2SurroundSCoordWeights;

spL=sum(SurroundLCoordWeights);
spM=sum(SurroundMCoordWeights);
spS=sum(SurroundSCoordWeights);

%Compute the DoGs
LDoG=CenterLCoordWeights-SurroundLCoordWeights;
MDoG=CenterMCoordWeights-SurroundMCoordWeights;
SDoG=CenterLCoordWeights-SurroundSCoordWeights;

%% SAVE TO OUTPUT

% GLOBAL
Cell.RetinaWeights=[pL pM pS];
Cell.CenterConeCount=[cnL cnM cnS];
Cell.SurroundConeCount=[snL snM snS];
Cell.CenterWeights=[cpL cpM cpS];
Cell.SurroundWeights=[spL spM spS];
Cell.CSStrengthRatio=CSStrengthRatio;

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
Cell.Center.CenterConeCoords=SortedTotalMosaic(1:ConesToCenter,:);

Cell.Center.LCones.Count=cnL;
Cell.Center.LCones.Coords=LCoords(1:ConesToCenter,:);
Cell.Center.LCones.EachWeight=CenterLCoordWeights;
Cell.Center.LCones.TotalWeight=cpL;

Cell.Center.MCones.Count=cnM;
Cell.Center.MCones.Coords=MCoords(1:ConesToCenter,:);
Cell.Center.MCones.EachWeight=CenterMCoordWeights;
Cell.Center.MCones.TotalWeight=cpM;

Cell.Center.SCones.Count=cnS;
Cell.Center.SCones.Coords=SCoords(1:ConesToCenter,:);
Cell.Center.SCones.EachWeight=CenterSCoordWeights;
Cell.Center.SCones.TotalWeight=cpS;

Cell.Center.CenterConeDistances=CenterX;
Cell.Center.CenterSigma=cSigma;
Cell.Center.CenterKernel=NormalizedCenterKernel;
Cell.Center.CenterKernelInField=fullfieldCenter;

% SURROUND
Cell.Surround.ConesToSurround=ConesToSurround;
Cell.Surround.SurroundConeCoords=SortedTotalMosaic(1:ConesToSurround,:);

Cell.Surround.LCones.Count=cnL;
Cell.Surround.LCones.Coords=LCoords(1:ConesToSurround,:);
Cell.Surround.LCones.EachWeight=SurroundLCoordWeights;
Cell.Surround.LCones.TotalWeight=cpL;

Cell.Surround.MCones.Count=cnM;
Cell.Surround.MCones.Coords=MCoords(1:ConesToSurround,:);
Cell.Surround.MCones.EachWeight=SurroundMCoordWeights;
Cell.Surround.MCones.TotalWeight=cpM;

Cell.Surround.SCones.Count=cnS;
Cell.Surround.SCones.Coords=SCoords(1:ConesToSurround,:);
Cell.Surround.SCones.EachWeight=SurroundSCoordWeights;
Cell.Surround.SCones.TotalWeight=cpS;

Cell.Surround.SurroundConeDistances=SurroundX;
Cell.Surround.SurroundSigma=cSigma;
Cell.Surround.SurroundKernel=NormalizedSurroundKernel;
Cell.Surround.SurroundKernelInField=fullfieldSurround;

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
