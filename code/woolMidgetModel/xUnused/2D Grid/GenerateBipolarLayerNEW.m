%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function Cell=GenerateBipolarLayer(Em)
%% GENERATE THE CONE MOSAIC
ConeMosaic=GenerateMosaic(Em);
%% GENERATE THE BIPOLAR MOSAIC
%Load a universal hex grid
cd('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid');
load('HexGrid');

%Jitter the cartesian coordinates to make a unique mosaic
pixel_coords(:,1)=jitter(pixel_coords(:,1));
pixel_coords(:,2)=jitter(pixel_coords(:,2));

%Determine the bipolar density given the eccentricity
BipolarsPerMM=ceil(((52170*exp(-1.057.*Em))+(15550*exp(-0.1605.*Em)))/(2.67)); %From 'cones per midget gc data' Excel file
DFRadius=-0.000565*(Em^2)+(0.01863*Em)+0.00146; %In mm. Polynomial fit of data from Fig 9, Watanabe and Rodieck JCN 1989

%Determine the aperture width needed to encompass the correct number of
%bipolars
r=1;
ConeCount=0;
while ConeCount<BipolarsPerMM
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
SortedTotalBipolars=sortrows([inc_pix_coordsADJ,distances],3);

%Determine each bipolar's nearest neighbors and use their average as an estimate of
%intercone distance across the retina
NearestNeighbors=zeros(length(SortedTotalBipolars),1);
HowManyNeighbors=1;
parfor q=1:length(SortedTotalBipolars)
    dist=sort(pdist2(SortedTotalBipolars(q,:),SortedTotalBipolars));
    NearestNeighbors(q)=mean(dist(2:HowManyNeighbors+1));
end
InterConeDistance=mean(NearestNeighbors);

%% BIPOLAR SUBSET
%Choose a central subset of the cone mosaic locations as your 
%substrate for the bipolar locations. This keeps processing speeds short...
BipolarSubstrate=SortedTotalBipolars(1:500,:);

%Identify the size of your ganglion cell and determine how many bipolars
%will fit within your DF
DFRadius=-0.000565*(Em^2)+(0.01863*Em)+0.00146;
BipolarToGanglionIndex=find(BipolarSubstrate(:,3)<DFRadius);
BipolarCells=BipolarSubstrate(find(BipolarSubstrate(:,3)<DFRadius),1:2);
NumberofBipolarInputs=length(BipolarCells);

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Location=BipolarCells(b,:);
end

%% CONE INPUTS TO EACH BIPOLAR

CSRadRatio=15;
CSStrengthRatio=.75;
ConesToCenter=7;
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Compute the distances of each cone from each bipolar and choose the ones
%that make up the inputs to that particular bipolar
MosaicL=~isnan(Cell.Global.LConeCoords(:,1));
MosaicM=~isnan(Cell.Global.MConeCoords(:,1));
ConeSubstrate=Cell.Global.AllConeCoords(1:5000,:);
LCones=MosaicL(1:length(ConeSubstrate));
MCones=MosaicM(1:length(ConeSubstrate));

for d=1:length(BipolarCells)
    BipolarDistances(:,d)=pdist2(ConeSubstrate(:,1:2),[BipolarCells(d,1),BipolarCells(d,2)]);
    [SortDistances(:,d), SortIndices(:,d)]=sort(BipolarDistances(:,d));
    BipolarConeCenterIndices(:,d)=SortIndices(1:ConesToCenter,d);
    BipolarConeSurroundIndices(:,d)=SortIndices(1:ConesToSurround,d);
end

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Total=ConeSubstrate(BipolarConeCenterIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,1)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,2)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.LWeight=sum(LCones(BipolarConeCenterIndices(:,b)))/length(LCones(BipolarConeCenterIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,1)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,2)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.MWeight=sum(MCones(BipolarConeCenterIndices(:,b)))/length(MCones(BipolarConeCenterIndices(:,b)));
end

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Total=ConeSubstrate(BipolarConeSurroundIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,2)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.LWeight=sum(LCones(BipolarConeSurroundIndices(:,b)))/length(LCones(BipolarConeSurroundIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,2)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.MWeight=sum(MCones(BipolarConeSurroundIndices(:,b)))/length(MCones(BipolarConeSurroundIndices(:,b)));
end

% %Assign a sigma to determine the shape of the kernels
% cSigma=max(CenterX)/3; %(we assume ~3 SDs make an RF radius)
% sSigma=max(SurroundX)/3;

for x=1:NumberofBipolarInputs
    CenterX=SortDistances(1:ConesToCenter,x);
    SurroundX=SortDistances(1:ConesToSurround,x);
    
    %Assign a sigma to determine the shape of the kernels
    cSigma=max(CenterX)/1; %(we assume ~1 SDs make an RF radius)
    sSigma=max(SurroundX)/1;

    %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
    %from the origin...
    CenterKernel(:,x)=exp((-1/(2*cSigma^2))*(abs(CenterX).^2));
    NormalizedCenterKernel(:,x)=CenterKernel/(sum(CenterKernel));
    SurroundKernel(:,x)=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
    NormalizedSurroundKernel(:,x)=(SurroundKernel/(sum(SurroundKernel)));
    
    %...and combine them to see the complete receptive field
    %DoG(:,x)=NormalizedCenterKernel-(CSStrengthRatio*NormalizedSurroundKernel);
end

PaddedNormalizedCenterKernel=padarray(NormalizedCenterKernel,size(NormalizedSurroundKernel,1)-size(NormalizedCenterKernel,1),'post');
DoG=PaddedNormalizedCenterKernel-(CSStrengthRatio*NormalizedSurroundKernel);

for h=1:NumberofBipolarInputs
    CenterLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalLWeight(:,h)=LCones(BipolarConeSurroundIndices(:,h)).*DoG(:,h);
    CenterMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*PaddedNormalizedCenterKernel(:,h);
    SurroundMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*(NormalizedSurroundKernel(:,h));
    TotalMWeight(:,h)=MCones(BipolarConeSurroundIndices(:,h)).*DoG(:,h);
end

%% WEIGHT BIPOLAR INPUTS TO THE PARASOL

%Apply Gaussian weight to each bipolar cell
BCdist=pdist2(BipolarCells,[0,0]);
BCsigma=max(BCdist)/1;
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

    

