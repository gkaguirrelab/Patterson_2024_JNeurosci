function ConeMosaic=GenerateMosaic_CrysS(Em)
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

%% SAVE TO OUTPUT

ConeMosaic.RetinaWeights=[pL pM pS];
ConeMosaic.ConesPerMM=ConeCount;
ConeMosaic.AllConeCoords=SortedTotalMosaic;
ConeMosaic.AllConeDistances=distances;
ConeMosaic.InterconeDistance=InterConeDistance;

ConeMosaic.LConeCoords=LCoords;
ConeMosaic.MConeCoords=MCoords;
ConeMosaic.SConeCoords=SCoords;