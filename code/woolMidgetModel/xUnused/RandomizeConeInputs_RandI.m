function [cnAll,snAll,cpAll,spAll]=RandomizeConeInputs_RandI(Em)
CellType=1;
%Em=5; %Eccentricity in mm
Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)


%Center cone count
SigmasToRFRadius=2; %Let's say that 2 SDs make up a RF radius
CSRadRatio=7; %Let's also say that a surround radius is ~7x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=1/2; %Finally, let's say that surround gain is ~1/2 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

DFradius=((0.00056*(Em)^2+0.0014*Em+0.0013)*1000)/223; %From Crook et al. Excel data. This is in degrees.
CenterSigma=DFradius/SigmasToRFRadius; %Here we assign a sigma for our cell given a particular eccentricity.
SurroundSigma=CSRadRatio*CenterSigma; %Ditto.

%L:M:S ratio
ConesPerMM=ceil(19890*Em^(-0.6331)); %Curcio 1990

pL=0.58;
pM=0.42;
pS=0;

LLIndex=0;
LUIndex=ceil(ConesPerMM*pL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+ceil(ConesPerMM*pM);
SLIndex=MUIndex+1;
SUIndex=ConesPerMM;

CenterAssign=randi(ConesPerMM,1,ConesToCenter);
cnL=length(find(CenterAssign<=LUIndex));
cnM=length(find(CenterAssign>MLIndex & CenterAssign<=MUIndex));
cnS=length(find(CenterAssign>SLIndex & CenterAssign<=SUIndex));

SurroundAssign=randi(ConesPerMM,1,ConesToSurround);
snL=length(find(SurroundAssign<=LUIndex));
snM=length(find(SurroundAssign>MLIndex & SurroundAssign<=MUIndex));
snS=length(find(SurroundAssign>SLIndex & SurroundAssign<=SUIndex));

% disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
% disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));
% disp(strcat(['Center weighting: ',num2str(cnL/(cnL+cnM+cnS)),' L, ',num2str(cnM/(cnL+cnM+cnS)),' M, ' ,num2str(cnS/(cnL+cnM+cnS)),' S']));
% disp(strcat(['Surround weighting: ',num2str(snL/(snL+snM+snS)),' L, ',num2str(snM/(snL+snM+snS)),' M, ' ,num2str(snS/(snL+snM+snS)),' S']));

cpL=cnL/(cnL+cnM+cnS);
cpM=cnM/(cnL+cnM+cnS);
cpS=cnS/(cnL+cnM+cnS);

spL=snL/(snL+snM+snS);
spM=snM/(snL+snM+snS);
spS=snS/(snL+snM+snS);

cnAll=[cnL cnM cnS];
snAll=[snL snM snS];

cpAll=[cpL cpM cpS];
spAll=[spL spM spS];