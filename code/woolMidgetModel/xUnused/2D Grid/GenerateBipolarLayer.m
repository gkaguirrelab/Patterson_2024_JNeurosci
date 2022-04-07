%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateBipolarLayer(Em)
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

%Compute the distances of each bipolar from the center (0,0) and sort the
%coordinates from closest to farthest
distances=pdist2(inc_pix_coordsADJ,[0,0]);
SortedTotalMosaic=sortrows([inc_pix_coordsADJ,distances],3);

%Determine each bipolar's nearest neighbors and use their average as an estimate of
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

%% BIPOLAR SUBSET
%Choose a central subset of the cone mosaic locations as your 
%substrate for the bipolar locations. This keeps processing speeds short...
BipolarSubstrate=SortedTotalMosaic(1:500,:);
MosaicLBipolar=MosaicL(1:length(BipolarSubstrate));
MosaicMBipolar=MosaicM(1:length(BipolarSubstrate));
NumberofBipolarInputs=7;

% %Adjust the density of bipolar cells (in terms of X intercone distances)
% BCDistance=InterConeDistance*2;
% 
% %Compute the distance of all bipolars to each bipolar
% NearestNeighbors=zeros(length(BipolarSubstrate));
% parfor q=1:length(BipolarSubstrate)
% dist=pdist2(BipolarSubstrate(q,:),BipolarSubstrate);
% NearestNeighbors(q,:)=dist;
% end
% 
% %Find the bipolars that are all at least BCDistance apart
% done=0;
% qualify=[];
% qualify(1)=1;
% vector=[];
% vector=find(NearestNeighbors(qualify(1):end,qualify(1))>BCDistance);
% qualify(2)=vector(1)+qualify(1)-1;
% 
% for x=2:length(BipolarSubstrate)
%     if x>length(qualify)
%         break
%     end
%     vector=find(NearestNeighbors(qualify(x):end,qualify(x))>BCDistance);
%     for y=1:length(vector)
%         if NearestNeighbors(vector(y)+qualify(x)-1,qualify(1:x-1))>BCDistance
%             qualify(x+1)=vector(y)+qualify(x)-1;
%             vector=[];
%             break
%         end
%     end
% end
% 
% %Log the bipolar coordinates and save only the first X, since this is the
% %assumed number of inputs to a parasol cell
% BipolarCells=BipolarSubstrate(qualify',1:2);
% BipolarCells=BipolarCells(1:NumberofBipolarInputs,:);

% For when bipolars are not 2x cones apart, but only 1
BipolarCells=BipolarSubstrate(1:NumberofBipolarInputs,1:2);

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Location=BipolarCells(b,:);
end

%% CONE INPUTS TO EACH BIPOLAR

CSRadRatio=6;
CSStrengthRatio=.75;
ConesToCenter=7;
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Compute the distances of each cone from each bipolar and choose the ones
%that make up the inputs to that particular bipolar
for d=1:length(BipolarCells)
    BipolarDistances(:,d)=pdist2(BipolarSubstrate(:,1:2),[BipolarCells(d,1),BipolarCells(d,2)]);
    [SortDistances(:,d), SortIndices(:,d)]=sort(BipolarDistances(:,d));
    BipolarConeCenterIndices(:,d)=SortIndices(1:ConesToCenter,d);
    BipolarConeSurroundIndices(:,d)=SortIndices(1:ConesToSurround,d);
end

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Total=BipolarSubstrate(BipolarConeCenterIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,1)=MosaicLBipolar(BipolarConeCenterIndices(:,b)).*(BipolarSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,2)=MosaicLBipolar(BipolarConeCenterIndices(:,b)).*(BipolarSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.LWeight=sum(MosaicLBipolar(BipolarConeCenterIndices(:,b)))/length(MosaicLBipolar(BipolarConeCenterIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,1)=MosaicMBipolar(BipolarConeCenterIndices(:,b)).*(BipolarSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,2)=MosaicMBipolar(BipolarConeCenterIndices(:,b)).*(BipolarSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.MWeight=sum(MosaicMBipolar(BipolarConeCenterIndices(:,b)))/length(MosaicMBipolar(BipolarConeCenterIndices(:,b)));
end

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Total=BipolarSubstrate(BipolarConeSurroundIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)=MosaicLBipolar(BipolarConeSurroundIndices(:,b)).*(BipolarSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,2)=MosaicLBipolar(BipolarConeSurroundIndices(:,b)).*(BipolarSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.LWeight=sum(MosaicLBipolar(BipolarConeSurroundIndices(:,b)))/length(MosaicLBipolar(BipolarConeSurroundIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)=MosaicMBipolar(BipolarConeSurroundIndices(:,b)).*(BipolarSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,2)=MosaicMBipolar(BipolarConeSurroundIndices(:,b)).*(BipolarSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.MWeight=sum(MosaicMBipolar(BipolarConeSurroundIndices(:,b)))/length(MosaicMBipolar(BipolarConeSurroundIndices(:,b)));
end

% %Assign a sigma to determine the shape of the kernels
% cSigma=max(CenterX)/3; %(we assume ~3 SDs make an RF radius)
% sSigma=max(SurroundX)/3;

for x=1:NumberofBipolarInputs
    CenterX=SortDistances(1:ConesToCenter,x);
    SurroundX=SortDistances(1:ConesToSurround,x);
    
    %Assign a sigma to determine the shape of the kernels
    cSigma=max(CenterX)/3; %(we assume ~3 SDs make an RF radius)
    sSigma=max(SurroundX)/3;

    %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
    %from the origin...
    CenterKernel(:,x)=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
    NormalizedCenterKernel(:,x)=CenterKernel/(sum(CenterKernel));
    SurroundKernel(:,x)=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
    NormalizedSurroundKernel(:,x)=(SurroundKernel/(sum(SurroundKernel)));
    
    %...and combine them to see the complete receptive field
    %DoG(:,x)=NormalizedCenterKernel-(CSStrengthRatio*NormalizedSurroundKernel);
end

PaddedNormalizedCenterKernel=padarray(NormalizedCenterKernel,length(NormalizedSurroundKernel)-length(NormalizedCenterKernel),'post');
DoG=PaddedNormalizedCenterKernel-(CSStrengthRatio*NormalizedSurroundKernel);

for h=1:NumberofBipolarInputs
    CenterLWeight(:,h)=MosaicLBipolar(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundLWeight(:,h)=MosaicLBipolar(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalLWeight(:,h)=MosaicLBipolar(BipolarConeSurroundIndices(:,h)).*DoG(:,h);
    CenterMWeight(:,h)=MosaicMBipolar(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundMWeight(:,h)=MosaicMBipolar(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalMWeight(:,h)=MosaicMBipolar(BipolarConeSurroundIndices(:,h)).*DoG(:,h);
end

%% WEIGHT BIPOLAR INPUTS TO THE PARASOL

%Apply Gaussian weight to each bipolar cell
BCdist=pdist2(BipolarCells,[0,0]);
BCsigma=max(BCdist)/3;
BCKernel=exp((-1/(2*BCsigma^2))*(abs(BCdist).^2));

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
spL=TotalLSurround/TotalSurround;
spM=TotalMSurround/TotalSurround;
cpS=0;
spS=0;
cpAll=[cpL cpM cpS];
spAll=[spL spM spS];

% Use these values for determining the center and surround sigmas for the
% DoG. They need to be converted to deg
MaxCenterDist=max(max(distances(1:ConesToCenter,:)));
MaxSurroundDist=max(max(distances(1:ConesToSurround,:)));
MaxCenterDistDeg=(MaxCenterDist*1000)/200;
MaxSurroundDistDeg=(MaxSurroundDist*1000)/200;

%% SAVE SOME SHIT

Cell.CenterWeights=[cpL cpM cpS];
Cell.SurroundWeights=[spL spM spS];
Cell.MaxCenterDistDeg=MaxCenterDistDeg;
Cell.MaxSurroundDistDeg=MaxSurroundDistDeg;

%% PLOT SOME SHIT

% figure;
% hold on;
% axis square;
% 
% plot(BipolarSubstrate(MosaicL(1:length(BipolarSubstrate)),1),BipolarSubstrate(MosaicL(1:length(BipolarSubstrate)),2),'ro');
% plot(BipolarSubstrate(MosaicM(1:length(BipolarSubstrate)),1),BipolarSubstrate(MosaicM(1:length(BipolarSubstrate)),2),'go');
% plot(BipolarCells(:,1),BipolarCells(:,2),'ko')
% 
% for p=1:NumberofBipolarInputs
%     plot(BipolarSubstrate(BipolarConeCenterIndices(:,p),1),BipolarSubstrate(BipolarConeCenterIndices(:,p),2),'wo')
%     plot(BipolarSubstrate(BipolarConeSurroundIndices(:,p),1),BipolarSubstrate(BipolarConeSurroundIndices(:,p),2),'wo')    
% end

    

