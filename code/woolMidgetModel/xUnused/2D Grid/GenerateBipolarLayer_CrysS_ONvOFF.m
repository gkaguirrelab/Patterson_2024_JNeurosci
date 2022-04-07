%function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=GenerateCell(Em)
function [ConeMosaic,BipolarLayer]=GenerateBipolarLayer_CrysS_ONvOFF(Em,CellType)
%% GENERATE THE CONE MOSAIC
ConeMosaic=GenerateMosaic_CrysS(Em);
%% GENERATE THE BIPOLAR MOSAIC
BipolarMosaic.AllBCCoords=ConeMosaic.AllConeCoords; %Midget bipolars are assumed to make 1:1 connections with cones at all eccentricities
BipolarMosaic.BCsPerMM=ConeMosaic.ConesPerMM;
BipolarMosaic.InterBCDistance=ConeMosaic.InterconeDistance;

%% BIPOLAR SUBSET
%Choose a central subset of the cone mosaic locations as your 
%substrate for the bipolar locations. This keeps processing speeds short...
BipolarCells=BipolarMosaic.AllBCCoords(1:500,1:2);
NumberofBipolarInputs=length(BipolarCells);

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Location=BipolarCells(b,:);
end

%% CONE INPUTS TO EACH BIPOLAR

CSRadRatio=10;
CSStrengthRatio=1;
ConesToCenter=1;
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

%Compute the distances of each cone from each bipolar and choose the ones
%that make up the inputs to that particular bipolar
MosaicL=~isnan(ConeMosaic.LConeCoords(:,1));
MosaicM=~isnan(ConeMosaic.MConeCoords(:,1));
MosaicS=~isnan(ConeMosaic.SConeCoords(:,1));

if length(ConeMosaic.AllConeCoords)>5000
    ConeSubstrate=ConeMosaic.AllConeCoords(1:5000,:);
else ConeSubstrate=ConeMosaic.AllConeCoords;
end

LCones=MosaicL(1:length(ConeSubstrate));
MCones=MosaicM(1:length(ConeSubstrate));
SCones=MosaicS(1:length(ConeSubstrate));

for d=1:length(BipolarCells)
    BipolarDistances(:,d)=pdist2(ConeSubstrate(:,1:2),[BipolarCells(d,1),BipolarCells(d,2)]);
    [SortDistances(:,d), SortIndices(:,d)]=sort(BipolarDistances(:,d));
    BipolarConeCenterIndices(:,d)=SortIndices(1:ConesToCenter,d);
    BipolarConeSurroundIndices(:,d)=SortIndices(1:ConesToSurround,d);
end

for b=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(b));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Total=ConeSubstrate(BipolarConeCenterIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.Distances=SortDistances(1:ConesToCenter,b);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,1)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.Locations(:,2)=LCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.LWeight=sum(LCones(BipolarConeCenterIndices(:,b)))/length(LCones(BipolarConeCenterIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,1)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.Locations(:,2)=MCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.MWeight=sum(MCones(BipolarConeCenterIndices(:,b)))/length(MCones(BipolarConeCenterIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.Locations(:,1)=SCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.Locations(:,2)=SCones(BipolarConeCenterIndices(:,b)).*(ConeSubstrate(BipolarConeCenterIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.SWeight=sum(SCones(BipolarConeCenterIndices(:,b)))/length(SCones(BipolarConeCenterIndices(:,b)));

    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Total=ConeSubstrate(BipolarConeSurroundIndices(:,b),1:2);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.Distances=SortDistances(1:ConesToSurround,b);
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,2)=LCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.LWeight=sum(LCones(BipolarConeSurroundIndices(:,b)))/length(LCones(BipolarConeSurroundIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,2)=MCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.MWeight=sum(MCones(BipolarConeSurroundIndices(:,b)))/length(MCones(BipolarConeSurroundIndices(:,b)));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.Locations(:,1)=SCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),1));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.Locations(:,2)=SCones(BipolarConeSurroundIndices(:,b)).*(ConeSubstrate(BipolarConeSurroundIndices(:,b),2));
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.SWeight=sum(SCones(BipolarConeSurroundIndices(:,b)))/length(SCones(BipolarConeSurroundIndices(:,b)));

end

% If the downstream cell is an ON cell, zero out the S-center bipolars
if CellType==1;
    for c=1:NumberofBipolarInputs
        CellLabel=strcat('Bipolar',num2str(c));
        BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.SWeight=0;
        BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.SWeight=0;
    end
end

% %Assign a sigma to determine the shape of the surround kernels. Remember
% the center is always 1 cone, so we can shape the ganglion gaussian later
% on

%Set up some parameters for the size and strength of the center and
%surround responses


for x=1:NumberofBipolarInputs
    CellLabel=strcat('Bipolar',num2str(x));
    SurroundX=SortDistances(1:ConesToSurround,x);
    
    %Assign a sigma to determine the shape of the kernels
    sSigma=max(SurroundX);

    %Apply a normalized rbfGaussian kernel to each cone as a function of its distance
    %from the bipolar center...
    SurroundKernel(:,x)=exp((-1/(2*sSigma^2))*(abs(SurroundX).^2));
    NormalizedSurroundKernel(:,x)=(SurroundKernel/(sum(SurroundKernel)));
    
    % H1 SURROUND COMPONENT
    H1StrengthRatio=.825; %H1 Cells make up ~82.5% of all horizontal synapses
    if BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.SWeight==1
        H1LMWeight=0; %If there are no L/M cones in the RF center, then there is no H1 contribution to the RF
    else
        H1LMWeight=1; %H1 cells only connect to L and M cones...
    end
    H1SWeight=0; % ...and never S cones :(
    
    % H2 SURROUND COMPONENT
    H2StrengthRatio=.175; %H2 cells make up 17.5% of all horizontal synapses
    H2LMWeight=1/6; %H2 signaling from a single L/M cone is 1/6 the strength of that from a single S cone
    H2SWeight=1;
    
    % Be sure to zero out the surrounds for any ON bipolars that *were*
    % S-center and are now nothing (because it's an ON cell) 
    if BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.LCones.LWeight+BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.MCones.MWeight+BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Center.SCones.SWeight==0
        % All the surrounds are = 0 because there is no center!
        H1LMWeight=0; %H1 cells only connect to L and M cones...
        H1SWeight=0; % ...and never S cones :(
        H2LMWeight=0; %H2 signaling from a single L/M cone is 1/6 the strength of that from a single S cone
        H2SWeight=0;
    end

    %For the H1 surround:
    H1SurroundLCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H1LMWeight*H1StrengthRatio;
    H1SurroundMCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H1LMWeight*H1StrengthRatio;
    H1SurroundSCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H1SWeight*H1StrengthRatio;

    %For the H2 surround:
    H2SurroundLCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H2LMWeight*H2StrengthRatio;
    H2SurroundMCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H2LMWeight*H2StrengthRatio;
    H2SurroundSCoordWeights=logical(BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.Locations(:,1)).*NormalizedSurroundKernel(:,x)*H2SWeight*H2StrengthRatio;

    %Total surround:
    SurroundLCoordWeights=H1SurroundLCoordWeights+H2SurroundLCoordWeights;
    SurroundMCoordWeights=H1SurroundMCoordWeights+H2SurroundMCoordWeights;
    SurroundSCoordWeights=H1SurroundSCoordWeights+H2SurroundSCoordWeights;
    
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.H1Surround=H1SurroundLCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.H1Surround=H1SurroundMCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.H1Surround=H1SurroundSCoordWeights;
    
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.H2Surround=H2SurroundLCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.H2Surround=H2SurroundMCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.H2Surround=H2SurroundSCoordWeights;
    
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.LCones.TotalSurround=SurroundLCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.MCones.TotalSurround=SurroundMCoordWeights;
    BipolarLayer.(matlab.lang.makeValidName(CellLabel)).Surround.SCones.TotalSurround=SurroundSCoordWeights;
    
    SurroundLWeight(:,x)=SurroundLCoordWeights;
    SurroundMWeight(:,x)=SurroundMCoordWeights;
    SurroundSWeight(:,x)=SurroundSCoordWeights;

end

BipolarLayer.SurroundArray.LConeWeights=SurroundLWeight;
BipolarLayer.SurroundArray.MConeWeights=SurroundMWeight;
BipolarLayer.SurroundArray.SConeWeights=SurroundSWeight;
BipolarLayer.SurroundArray.Indices=BipolarConeSurroundIndices;
